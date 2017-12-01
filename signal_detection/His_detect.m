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
%data = s.data;
data0 = s.PATIENT1SINUSRHYTHMNUMBERSONLY;
%data0 = s.PATIENT1SUPRAVENTRICULARTACHYCARDIANumbersOnly;

%Define the amount of data that we want to use for parameter learning and
%detection.
begin_time = 0.0;
end_time = 20; %second
data = data0(:,16);
%{
%filter to remove DC bias
b = fir1(1000,2.5/Fs,'high');
b2 = fir1(1000,150/Fs);
% use only the HRA channel from the data
data = filter(b2,1,filter(b,1,data(:,16)));
data = data(begin_time*Fs+1:end_time*Fs+1,:);
%}
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
%%
his = data0(:,10:14);
plot(his(:,1));
plot(his(:,2));
plot(his(:,3));
plot(his(:,4));
plot(his(:,5));

temphis = his(aind(23):aind(24),:);
figure();
plot(temphis(:,5))
%% negative peak finder
step = 10;
vwidth= 100;
HisInfo = zeros(2,5);
finalHis = cell(1,5);
vstamps = zeros(1,5);
for i = 1:5
    optimal_vstamp = length(temphis(:,i));
    Npeakval = -7000;
    [~,vindex] = max(temphis(:,i));
    for j = 1:600
        nNpeak = sum(temphis(:,i) <= Npeakval);
        nNindex = find(temphis(:,i) <= Npeakval);
        if (vindex-min(nNindex) <= vwidth) & (nNpeak <= 40)
            optimal_vstamp = min(nNindex);
            optimal_cutoff = Npeakval;
        end
        Npeakval = Npeakval+step;
        
    end
    vstamps(i) = optimal_vstamp;
    figure();
    finalHis(1,i) = {temphis(1:optimal_vstamp,i).'};
    hold on;
    plot(temphis(:,i));
    plot(temphis(1:optimal_vstamp,i).'); 
    title(['Final his chunk',num2str(i)])
    hline = refline([0 optimal_cutoff]);
    hline.Color = 'r';
    
    %% 
    %hold off;
    %figure();
    %finalHis = temphis(1:optimal_vstamp,i);
    %plot(finalHis);
    %}
end

%% his finder

vind = min(vstamps);
HisWSize = 10;
for i = 1:5
    Hisdiff = padarray(diff(temphis(1:vind-20,i).'),[0 1],'pre');
    HisFluc = abs(((temphis(1:vind-20,i).^2).').*Hisdiff);
    HisWSum = zeros(1,length(HisFluc)-HisWSize+1);
    for j = 1:length(HisWSum)
        HisWindow = HisFluc(j:min(length(HisFluc),j+HisWSize-1));
        HisWSum(j) = sum(HisWindow);
    end
    [V,I] = max(HisWSum);
    HisInfo(:,i) = [V,floor(I+HisWSize/2)];

end
figure();
hold on;

plot(temphis(:,1));
plot(temphis(:,2));
plot(temphis(:,3));
plot(temphis(:,4));
plot(temphis(:,5));
legend('1','2','3','4','5')
[HisV,HisI] = max(HisInfo(1,:));
hisy = cell2mat(finalHis(1,HisI));
plot(HisInfo(2,HisI),hisy(HisInfo(2,HisI)),'b*');

%}



