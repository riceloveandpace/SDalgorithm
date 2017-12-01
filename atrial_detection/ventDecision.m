function ventDecision

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

close all
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
data = data(:,10);
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
[ds.thresh,ds.flip,ds.length]=atrialParamLearning(datalearn);

ds.beatDelay = 0; %Tracks amount of time since last ventricular beat.
ds.beatFallDelay = 0;%Tracks amount of time since last falling edge of ventricular beat.
ds.PostVARP = 250;%Minimum time between ventricular then atrial beat.
ds.PreVARP = 20; % check!!!!!
ds.PostAVRP = 100;%Minimum time between atrial then ventricular beat.
ds.PreAVRP = 20; % check
ds.VV = 350;

%these are used for doing real time detection
ds.recentBools = zeros(1,ds.length); %Binary value indicating whether recent values have exceeded v_thresh
ds.last_sample_is_sig = false;%Flag indicating whether last sample was a V beat
%
%These variables are just used for testing and visualization, not actually
%used in the algorithm.
ds.PeakInd = [];
ds.recentdatapoints = zeros(1,ds.VV);

%This loop models real time data acquisition in an actual hardware system.
for i = 1:numSamples
    
    %Increment time since last atrial beat d.
    ds.recentdatapoints = [ds.recentdatapoints(2:end) data(i)];
    %get next datapoint and add to buffer
    
    %increment all delays when considering each sample in real time.
    ds.beatDelay = ds.beatDelay + 1;
    ds.beatFallDelay = ds.beatFallDelay + 1;
    
    %Perform beat detection with the knowledge of the new sample.
    ds = singlePeakFinder(i,ds);
end

figure
d = ds;
hold on;
h=plot(data,'k');
a=plot([0 d.PeakInd], d.thresh*d.flip, 'or'); h=[h a(1)];
title(['Channel ' num2str(i)],'Fontsize',14)
legend(h,{'data','ventricular peaks','atrial peaks'},'Fontsize',14)
end
