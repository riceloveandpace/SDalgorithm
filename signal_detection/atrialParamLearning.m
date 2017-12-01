function [ height, flip, a_length] = atrialParamLearning(data)

    %learn the appropriate lengths of atria, as well as
    %which should be detected first, for + and - data
    lengths = LearnLengthsnew(data);
    
    visualize_bs = true; % Visualize the function that LearnHeightRange uses and the algorithm's performance
    flip = [+1, -1]; % signs corresponding to the two flip indices
    
    %The ranges of magnitudes returned by LearnHeightRange for each chamber and sign
    cutoffs = zeros(2,2); % Dimensions: flip, start-end
    
    %Where in the middle of each range we should pick the magnitude
    %threshold; each row should sum to 1
    mid_finders = [0.5 0.5; 0.5 0.5]; % Dimensions: chamber, start-end
    
    if visualize_bs
        figure
        hold on
    end
    
    % + and - data
    for f = 1:2
        flip(f)
        %Learn the cutoffs for whichever chamber we should do first for this sign of the data
        cutoffs(f,:) = LearnHeightRange(flip(f)*data, lengths(f));
        
        %Find peak times and zero out the region around them
        ndata = data;
        t_blank = 100; l_blank = 100; %how many samples to zero out on either side
        %pick exact threshold from the middle of the range according to mid_finders
        thresh = sum(squeeze(cutoffs(f,:)).*mid_finders);

        [~, rising_edges, falling_edges] = CountPeaks(flip(f)*data > thresh(f), lengths(f)); %find peak times
        %zero out surrounding region
        for i = 1:length(rising_edges)
            ndata(max(1,(rising_edges(i)-t_blank)):min(end,(falling_edges(i)+l_blank))) = 0;
        end
        
        if visualize_bs
                color = 'r';
                peak = 'atrial';
            %end
            
            xs = linspace(0, max(data), 200);
            ys=zeros(200,1);
            for i = 1:200
                [beats, ~, ~] = CountPeaks(flip(f)*data > xs(i), lengths(f));
                ys(i) = beats;
            end
            h(1) = plot(flip(f)*xs, ys,color);
            
            % axis([0 3000 0 200])
            xlabel('Threshold (mV)','FontSize',18)
            ylabel('Number of Beats','FontSize',18)
            title('Beats(threshold)','Fontsize',18)
            xs = linspace(0, max(flip(f)*ndata), 200);
            ys= zeros(200,1);
            for i = 1:200
                [beats, ~, ~] = CountPeaks(flip(f)*ndata > xs(i), lengths(f));
                ys(i) = beats;
            end
            h(2) = plot(flip(f)*xs, ys, color);
        end
    end
    
    %pick which signs to use according to which one returned a wider range in LearnHeightRange
    [~,af] = max(cutoffs(:,2)-cutoffs(:,1));
    %pick exact thresholds from the middle of the ranges according to mid_finders
    height = sum(cutoffs(af,:).*mid_finders(af));
    a_length = lengths(af);
    flip = flip(af);
end


% Learn an appropriate threshold for beat detection 
% by searching (the function that maps threshold to # of beats detected in data)
% for a maximal-length flat segment. This corresponds to a large range of thresholds 
% that give the same number of beats, and we've found that those thresholds are good.
function [flatcutoffs] = LearnHeightRange(data, minlength)
    if (minlength==0) % length learning failed, so threshold learning doesn't make sense
        flatcutoffs = [0 0];
        return
    end
    npts = 20; % How many samples to take in each iteration of our recursive search
    sample_rate = 1000;
    minbeats = length(data) / sample_rate * 10 / 60; %10 bpm
    maxbeats = length(data) / sample_rate * 150 / 60; %200bpm

    max_th = max(data); % upper bound for this iteration of the search
	% but try two other uppper bounds. Having a better upper bound makes the search better
	if CountPeaks(data > max(data) / 4, minlength) < minbeats
		max_th = max(data)/4 % We already know this is too big, so no need to check anything larger
	elseif CountPeaks(data > max(data) / 2, minlength) < minbeats
		max_th = max(data)/2 
	end
    
    ths = linspace(0, max_th, npts); % thresholds to try in this iteration
    beats = zeros(1,npts);
    
   
    for i=1:3 % refinement round
        beats = zeros(1,npts);
        for j = 1:length(ths)
            % Potential optimization: We can assume the function is monotonic, so if beats(j) =
            % beats(k), we don't need to compute anything else between j
            % and k.
            beats(j) = CountPeaks(data > ths(j), minlength);
			% as j increases, ths(j) increases, so beats(j) decreases. 
			% This means if beats(j) < minbeats, then beats(j+1) < minbeats, etc. 
            if beats(j)<minbeats
				% don't touch the rest. They'll be 0, which we treat the same as any value < minbeats
                break
            end
        end
        last_valid = find(beats > minbeats, 1, 'last');
        first_valid = find(beats < maxbeats, 1, 'first');
		% Could we find any thresholds that were reasonable?
        if ( isempty(last_valid) || isempty(first_valid) ||(last_valid<=first_valid))
            flatcutoffs = [0 0];
            return
        end
        
        beats = beats(first_valid:last_valid); % just take the valid ones
        ths = ths(first_valid:last_valid);
        derivs = abs(beats(2:end)-beats(1:(end-1))); % consecutive differences
        %derivs = abs(beats(2:last_valid)-beats(1:(last_valid-1)));
        [minder, idx] = min(derivs); % which segment (between points we sampled) is the flattest?
        if minder == 0
            % was it actually flat? If so, find the longest strech of zeros
            flats = derivs == 0;
            iidx = -1;
            currentBest = -1;
            cbStart = -1;
            cbEnd = -1;
            for j = 1:length(derivs)
                if flats(j) && iidx  == -1 % extending a flat?
                    iidx = j;
                end
                if flats(j) && j - iidx > currentBest
                    currentBest = j - iidx;
                    cbStart = iidx;
                    cbEnd = j;
                end
                if ~flats(j)
                    iidx = -1;
                end 
            end
            idx = cbStart;
            endidx = cbEnd+1;
        else
			% otherwise, take the flattest segment
            endidx = idx+1;
        end
		% compute the thresholds for the next round
        nmin = min(ths(idx),ths(endidx));
        nmax = max(ths(idx),ths(endidx));
        ths = linspace(nmin, nmax, npts);
        
        beats1 = beats(idx);
        beats2 = beats(endidx);
		% if we have already reached a flat interval, that's good enough. We're done.
        if beats1 == beats2
            break
        end
		% we know these values, so we actually don't need to recompute them.
        % beats(1) = beats1;
        % beats(end) = beats2;
    end
    % Technically, this might not be the full extent of the interval because of how we sampled.
	% We could try to expand it, but that's not particularly necessary as far as I can tell.
    flatcutoffs = [ths(1), ths(end)];
end

function [lengths] = LearnLengthsnew(data)

%featurize the data into a list of peaks with [width, height] values by
%locating edges at the pseudo-steepest points in the waveform (extrema in
%the weigthed derivative), and then featurizing peaks as appropriate pairs of
%rising and falling edges

extreme_over = 100; 
ddatadt = [0; diff(data)].^2.*sign([0; diff(data)]).*data.*sign(data);
ddatadt = ddatadt * max(data)/max(ddatadt);
lengths = [];
for flip = [+1 -1]
    wall_times = [];
    wall_steeps = [];
    for i = 1:length(ddatadt)
        nearby = ddatadt(max(1,i-extreme_over):min(end,i+extreme_over)).*(data(max(1,i-extreme_over):min(end,i+extreme_over))*flip>0);
        if flip*data(i)>0 && (all(ddatadt(i) >= nearby) || all(ddatadt(i) <= nearby))
            wall_times = [wall_times i];
            wall_steeps = [wall_steeps ddatadt(i)];
        end
    end
    
    if (false)
        figure
        plot(data)
        hold on
        plot(ddatadt)
        stem(wall_times, wall_steeps)
        xlabel('time (samples)','Fontsize',14)
        legend({'input waveform', 'weighted derivative of input', 'steepest points'},'Fontsize',14)
        title('Finding the steepest points','Fontsize',18)
    end

    peak_lengths = [];
    peak_steeps = [];
    peak_heights = [];
    peak_highests = [];
    for i = 2:length(wall_times)
        len = wall_times(i)-wall_times(i-1);
        highest = sign(data(wall_times(i-1)))*max(abs(data(wall_times(i-1):wall_times(i))));
        height = highest - (data(wall_times(i-1))+data(wall_times(i)))/2;
        if ( wall_steeps(i-1)*wall_steeps(i) < 0 && flip*data(wall_times(i-1))>0 && flip*data(wall_times(i))>0 && flip*wall_steeps(i-1)>0 && flip*wall_steeps(i)< 0 && len < 75)% && s*height > abs(median(data)) )
            peak_lengths = [peak_lengths; len];
            peak_steeps = [peak_steeps; wall_steeps(i-1) - wall_steeps(i)];
            peak_heights = [peak_heights; height];
            peak_highests = [peak_highests; highest];
        end
    end
    
    %remove outliers
	% simpler version:
    % outliers = abs(zscore([peak_lengths peak_heights peak_steeps]))>3;
    zscores = zscore([peak_lengths peak_heights peak_steeps]);
    zscorenorms = sqrt(sum(zscores.^2,2));
    outliers = zscorenorms > 2.5;
    to_removes = [];
    for i=1:length(peak_lengths)
        if any(outliers(i,:))
            to_removes = [to_removes i];
        end
    end
    if (length(to_removes)/length(peak_lengths)>1/10)
        disp('Excluding a lot of outliers')
    end
    peak_lengths(to_removes) = [];
    peak_steeps(to_removes) = [];
    peak_heights(to_removes) = [];
    peak_highests(to_removes) = [];
    
    if (length(peak_lengths)<=10)
        lengths = [lengths 0];
        continue
    end
    
    lengths = [lengths round(mean(peak_lengths))];
end
lengths
end
