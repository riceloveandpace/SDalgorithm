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
addpath('test_data')
addpath('test_data/sept29_2016_test/')
addpath('test_data/PhisioBank_iaf/')
%addpath('test_data/FebruaryData/')
% each .mat file contains a struct for data from a test
% currently, each struct has two fields:
% sampling_rate: the sampling rate of the data in Hz
% data: a NxC matrix where each column represents a channel and each row
%       represents a time point

s = load('ep1SR.mat');

%Set the sampling rate and data that we will be working with.
%Fs = s.Fs;
Fs = 1000;
data = s.PATIENT1SINUSRHYTHMNUMBERSONLY;
if in1 == 'ven1'
    data = data(:,17);
elseif in1 == 'atr1'
    data = data(:,16);
elseif in1 == 'atr2'
    data = data(:,15);
elseif in1 == 'ven2'
    data = data(:,18);
end
%size(data)
%Define the amount of data that we want to use for parameter learning and
%detection.
begin_time = 0.0;
end_time = 7; %second
datalearn = data(begin_time*Fs+1:end_time*Fs+1,:);
%data = [data(:,1) data(:,2) data(:,10:14)];
[numSamples, numChannels] = size(data);
size(data)
%~~~~~~~~~~~~~~~~~~~~~
%Channel Parameters
%~~~~~~~~~~~~~~~~~~~~~
ds = struct();

%Compute all of the detection parameters for this channel using
%LearnParameters (LearnParameters will also call LearnLengths.
[ds.thresh,ds.flip,ds.len]=atrialParamLearning(datalearn);
% for simplicity I referred to the non-real-time algorithm
[aindlearn] = atrial_peak_finder(ds, datalearn); 
% learn the energy threshold
j = 1; ennoise = [];
windlen = 8;
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
noiseavg = mean(ennoise);
ds.noiseavg = noiseavg;

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
%%
%This loop models real time data acquisition in an actual hardware system.
storlen = 120;
ds.storen = zeros(1,storlen); % this need to be stored for previous samples !!!
ds.startind = [];
ds.endind = [];
ds.findEnd = 'f';
enStor = 't';
ennew = []; enallnew = [];
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
    % Perform start detection: need some memory
    if ds.last_sample_is_sig % then go and find start
        ds.findEnd = 't';
        tempen = ds.storen(end);
        k = 0;
        ds.storen;
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

figure
d = ds;
hold on;
h=plot(data,'k');
a=plot([0 d.PeakInd], d.thresh*d.flip, 'or'); h=[h a(1)];
plot([0 ds.startind'],500,'xr')
plot([0 ds.endind'],500,'xb')
en = [];
for i = 1:length(data)-5
    tempdat = data(i:i+5);
    meannoise = mean(tempdat);
    en = [en; sumabs(tempdat)-5*abs(meannoise)];
end
plot(en)
plot([0,60000],[noiseavg,noiseavg])
%%
if in1 == 'ven'
    title('Ventricular Channel','Fontsize',14)
elseif in1 == 'atr'
    title('Atrial Channel','Fontsize',14)
legend(h,{'data','peaks'},'Fontsize',14)
end
plot([0,atrall(:,1)'],400,'or')
plot([0,atrall(:,2)'],400,'ob')
