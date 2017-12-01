% atrial detection: this is a simplified version of multisite_detect and
% functions on the HRA channel to detect the atrial signal
% not real-time
%Input:
%s - Structure that will contain the data, with multiple channels, and
%appropriate sampling rate.
%Output:
%Plots that will showcase detection of ventricular and atrial beats using
%parameters learned for each channel.
%close all
clear
%Paths for data that we have been using for the sake of testing.
addpath('test_data/sept29_2016_test/');
addpath('test_data');
%s = load('ep1SVT.mat');
s = load('ep1SR.mat');
%s = load('ep1TerminatingTachy.mat');
%Set the sampling rate and data that we will be working with.
%Fs = s.Fs; %
Fs = 1000;
% data = s.data;
data = s.PATIENT1SINUSRHYTHMNUMBERSONLY;

%Define the amount of data that we want to use for parameter learning and
%detection.
begin_time = 0.0;
end_time = 7; %second
% use only the HRA channel from the data
data = data(:,16);
%data = filter(b2,1,filter(b,1,data(:,16)));
datalearn = data(begin_time*Fs+1:end_time*Fs+1,:);
datadown = datalearn(1:5:end);
[numSamples, numChannels] = size(data);

ainds = zeros(numChannels, 1);
%Compute all of the detection parameters for this channel using LearnParameters
[d(1).thresh, d(1).flip,d(1).alength] = atrialParamLearning(datalearn(:,1));
[aindlearn] = atrial_peak_finder(d(1), datalearn(:,1));
% learn the energy threshold
j = 1; ennoise = [];
for i = 1:length(datalearn(:,1))-d(1).alength % (+-80) out of the atrial peak count as noise
    if j <= length(aindlearn)
        if i < aindlearn(j)-80 % i.e. smaller than a peak
            ennoise = [ennoise;sumabs(data(i:i+d(1).alength))];
        elseif i > aindlearn(j) - 80 && i < aindlearn(j) + 80
            disp('in a beat')
        elseif i == aindlearn(j) + 80
            ennoise = [ennoise;sumabs(data(i:i+d(1).alength))];
            j = j + 1;
        end
    else
        ennoise = [ennoise;sumabs(data(i:i+d(1).alength))];
    end
end
noiseavg = mean(ennoise);

%Beat detection w/ the learnt parameters
[aind1] = atrial_peak_finder(d(1), data(:,1));
figure; hold on;
plot(data,'b');
plot(padarray(aind1.',[0 1],'pre'), d.thresh*d.flip, 'xg');
title('Channel 16','Fontsize',18)
%legend({'input waveform','detected ventricles','detected atria'},'Fontsize',18)
xlabel('time (samples)','Fontsize',14)

% start/stop detection
j = 1; en = [];
for i = 1:length(data)-d(1).alength
    temppeak = aind1(j);
    tempdat = data(i:i+d(1).alength);
    en = [en; sumabs(tempdat)];
end
startind = []; endind = [];
for j = 1:length(aind1)
    for i = 1:length(data)-d(1).alength
        if i == aind1(j) % if the current timestamp is a peak
            tempen = en(i);
            k = 0;
            while tempen >= noiseavg
                % search back for start
                k = k - 1;
                tempen = en(i+k);
                if tempen < noiseavg
                    startind = [startind;i+k];
                end
            end
            tempen = en(i); l = 0;
            while tempen >= noiseavg
                % search forward for end
                l = l + 1;
                tempen = en(i+l);
                if tempen < noiseavg
                    endind = [endind;i+l];
                end
            end
        end
    end
end
plot(startind,0, 'xr','MarkerSize',10);
plot(endind,0,'xk','MarkerSize',10);
