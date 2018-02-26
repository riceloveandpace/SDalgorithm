%%function ventDecision

%This script simulates our algorithm in real time.  It performs parameter
%learning on each channel of data.  In addition, it will perform beat
%detection on each of the channels.
%
%Input:
%s - Structure that will contain the data, with multiple channels, and
%appropriate sampling rate.
%
%Output:
%Plots that will showcase ventricular beat detection.
in1 = input('Enter channel: "atr", or "ven"','s');

%close all
addpath('PacingData')
addpath('EPStudys')

% each .mat file contains a struct for data from a test
% currently, each struct has two fields:
% sampling_rate: the sampling rate of the data in Hz
% data: a NxC matrix where each column represents a channel and each row
%       represents a time point

%s = load('IData1000Hz.mat');
%data = s.IData1000Hz;

%s = load('ep1SR.mat');
%data = s.PATIENT1SINUSRHYTHMNUMBERSONLY;

s = load('ep1SVT.mat');
data = s.PATIENT1SUPRAVENTRICULARTACHYCARDIANumbersOnly;

%Set the sampling rate and data that we will be working with.
%Fs = s.Fs;
Fs = 1000;
if in1 == 'ven1'
    data = data(:,14); %6
elseif in1 == 'atr1'
    data = data(:,15); %1
elseif in1 == 'atr2'
    data = data(:,16);
elseif in1 == 'ven2'
    data = data(:,18);

end
%size(data)
%Define the amount of data that we want to use for parameter learning and
%detection.
begin_time = 0.0;
end_time1 = 3; %second
end_time2 = 10;
end_time3 = 20;
datalearn1 = data(begin_time*Fs+1:end_time1*Fs+1,:);
datalearn2 = data(end_time1*Fs+1:end_time2*Fs+1,:);
datalearn3 = data(end_time2*Fs+1:end_time3*Fs+1,:);
%data = [data(:,1) data(:,2) data(:,10:14)];
[numSamples, numChannels] = size(data);
size(data);
%~~~~~~~~~~~~~~~~~~~~~
%Channel Parameters
%~~~~~~~~~~~~~~~~~~~~~
ds = struct();
ds.PeakInd = [];
ds.lastPI = 0;

ds.VV = 400;
%ds.maxennoise = 0;
%Compute all of the detection parameters for this channel using
%LearnParameters (LearnParameters will also call LearnLengths.
[ds]=atrialParamLearning(datalearn1,datalearn2,datalearn3,ds);
% for simplicity I referred to the non-real-time algorithm
%[aindlearn] = atrial_peak_finder(ds, datalearn); 
% learn the energy threshold
%j = 1; ennoise = [];
%windlen = 5;
%ds.beatDelay = 0; %Tracks amount of time since last ventricular beat.
%ds.beatFallDelay = 0;%Tracks amount of time since last falling edge of ventricular beat.

%these are used for doing real time detection
%ds.recentBools = zeros(1,ds.len); %Binary value indicating whether recent values have exceeded v_thresh
%ds.last_sample_is_sig = false;%Flag indicating whether last sample was a V beat
%
%These variables are just used for testing and visualization, not actually
%used in the algorithm.
%{
for i = 1:length(datalearn(:,1))-windlen%ds.len % (+-80) out of the atrial peak count as noise
    if j <= length(aindlearn)
        if i < aindlearn(j)-80 % i.e. smaller than a peak
        %    if i > aindlearn(j) - 500 % start
            ennoise = [ennoise;sumabs(data(i:i+windlen))-windlen*abs(mean(data(i:i+windlen)))];%ds.len))];
        %    else % end
        %    ennoise2 = [ennoise2;sumabs(data(i:i+windlen))-windlen*abs(mean(data(i:i+windlen)))];%ds.len))];
        %    end
        elseif i > aindlearn(j) - 80 && i < aindlearn(j) + 80
          %  disp('in a beat')
        elseif i == aindlearn(j) + 80
            ennoise = [ennoise;sumabs(data(i:i+windlen))-windlen*abs(mean(data(i:i+windlen)))];%ds.len))];
            j = j + 1;
        end
    else
        ennoise = [ennoise;sumabs(data(i:i+windlen))-windlen*abs(mean(data(i:i+windlen)))];%ds.len))];
    end
end
noiseavg = max(ennoise);
ds.noiseavg = noiseavg;

%{

ds.beatDelay = 0; %Tracks amount of time since last ventricular beat.
ds.beatFallDelay = 0;%Tracks amount of time since last falling edge of ventricular beat.
ds.VV = 350;

%these are used for doing real time detection
ds.recentBools = zeros(1,ds.len); %Binary value indicating whether recent values have exceeded v_thresh
ds.last_sample_is_sig = false;%Flag indicating whether last sample was a V beat
%
%These variables are just used for testing and visualization, not actually
%used in the algorithm.
ds.PeakInd = [];
ds.recentdatapoints = zeros(1,ds.VV);
%}

%}
%%
%This loop models real time data acquisition in an actual hardware system.
%storlen = 120;
%ds.storen = zeros(1,storlen); % this need to be stored for previous samples !!!
ds.startind = [];
ds.endind = [];
ds.findEnd = 'f';
ds.findPeak = 'f';
%enStor = 't';
ds.lastPI = 0;
%ennoise = 0;
%ennew = []; enallnew = [];

%detection of peaks based on learned parameters  && (i > ds.lastPI + ds.VV - 80)

for i = 1:numSamples
    
    k = data(i);
    if (abs(k) < ds.thresh) && (ds.findEnd == 'f') && (ds.findPeak == 'f')
        
        if abs(k) > ds.noiselvl%ds.noiselvl %found start index, ready to detect end point
            
            ds.startind = [ds.startind; i];
            ds.findPeak = 't';
        end
        
    elseif (k > ds.thresh) && (ds.findPeak == 't') && (ds.findEnd == 'f')
            %found peak, ready to detect peak index
            ds.PeakInd = [ds.PeakInd; i];
            ds.lastPI = i;
            ds.findEnd = 't';
            ds.findPeak = 'f';
            
            %found end, start over for a new beat
    elseif (abs(k) < ds.noiselvl) && (i > ds.lastPI + 80) && (ds.findEnd == 't') && (ds.findPeak == 'f')
                %found start index, ready to detect end point
                ds.endind = [ds.endind; i];
                ds.findPeak = 'f';
                ds.findEnd = 'f';
               
        
   
        
        
    end
       
    
end



%{
for i = 1:numSamples-5
    %Increment time since last atrial beat d.
    ds.recentdatapoints = [ds.recentdatapoints(2:end) data(i)];
    
    % store the energy every 2 timestamps
   % if enStor == 't'
   tempmean = mean(ds.recentdatapoints(ds.VV-windlen:ds.VV));
   ennew = sumabs(ds.recentdatapoints(ds.VV-windlen:ds.VV))-windlen*abs(tempmean);%ds.len+1:ds.VV))
   ds.storen = [ds.storen(2:end) ennew];
   %     enStor = 'f';
  %  else
  %      enStor = 't';
  %  end
    %increment all delays when considering each sample in real time.
    ds.beatDelay = ds.beatDelay + 1;
    ds.beatFallDelay = ds.beatFallDelay + 1;
    
    
    %Perform beat detection with the knowledge of the new sample.
    ds = singlePeakFinder(i,ds);
    
    %Found peak index, can do pacing decisions
    
    
    
    
    % Perform start detection: need some memory
    if ds.last_sample_is_sig % then go and find start
        ds.findEnd = 't';
        tempen = ds.storen(end);
        k = 0;
        while tempen >= noiseavg
            % search back for start
            k = k - 1;
            tempen = ds.storen(storlen+k);
            
            if tempen < noiseavg
                ds.startind = [ds.startind;i+k-windlen];
            end
        end
    end
    % wait to perform end detection
    if ds.findEnd == 't' % for each subsequent datapoint
        if ds.storen(end) < noiseavg
            ds.endind = [ds.endind;i+windlen];
            ds.findEnd = 'f';
        end
    end
    
end

ds.startend = [ds.startind ds.endind];
%}
%figure
figure(1);
%d = ds;
hold on;
plot(data,'k');
plot( ds.PeakInd, ds.thresh, 'or');% h=[h a(1)];
plot( ds.startind,500,'xr')
plot( ds.endind,500,'xb')
%detectStartend = ds.startend;
xlabel('timestamp')
ylabel('Magnitude')
%save('detect_ep1SR_ventall_Ch18.mat','detectStartend')
if in1 == 'ven1'
     title('Ventricular Channel 1','Fontsize',14)

elseif (in1 == 'atr1' ) 
     title('Atrial Channel','Fontsize',14)
elseif (in1 == 'atr2')
     title('Atrial Channel 2','Fontsize',14)
end

% end
%plot([0,atrall(:,1)'],400,'or')
%plot([0,atrall(:,2)'],400,'ob')
%}
