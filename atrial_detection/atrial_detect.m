% atrial detection: this is a simplified version of multisite_detect and
% functions on the HRA channel to detect the atrial signal
% not real-time
%Input:
%s - Structure that will contain the data, with multiple channels, and
%appropriate sampling rate.
%Output:
%Plots that will showcase detection of ventricular and atrial beats using
%parameters learned for each channel.
close all
clear
%Paths for data that we have been using for the sake of testing.
addpath('test_data/sept29_2016_test/');
addpath('test_data');
%s = load('ep1SVT.mat');
s = load('ep1SR.mat');
%Set the sampling rate and data that we will be working with.
%Fs = s.Fs; %
Fs = 1000;
% data = s.data;
data = s.PATIENT1SINUSRHYTHMNUMBERSONLY;
%data = s.PATIENT1SUPRAVENTRICULARTACHYCARDIANumbersOnly;

%Define the amount of data that we want to use for parameter learning and
%detection.
begin_time = 0.0;
end_time = 20; %second
%filter to remove DC bias
b = fir1(1000,2.5/Fs,'high');
b2 = fir1(1000,150/Fs);
% use only the HRA channel from the data
data = filter(b2,1,filter(b,1,data(:,16)));
data = data(begin_time*Fs+1:end_time*Fs+1,:);
[numSamples, numChannels] = size(data);

ainds = zeros(numChannels, 1);
%Compute all of the detection parameters for this channel using LearnParameters
[d(1).thresh, d(1).flip,d(1).alength] = atrialParamLearning(data(:,1));

%Beat detection w/ the learnt parameters
aind = atrial_peak_finder(d(1), data(:,1));
figure; hold on;
plot(data,'b');
plot(padarray(aind.',[0 1],'pre'), d.thresh*d.flip, 'xk');
title('Channel 16','Fontsize',18)
%legend({'input waveform','detected ventricles','detected atria'},'Fontsize',18)
xlabel('time (samples)','Fontsize',14)
