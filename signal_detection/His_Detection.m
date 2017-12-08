
aindlearn = ds.PeakInd;
HisInfo = zeros(length(aindlearn),2);
dh = struct();
dh.NumberOfChannel = 5;
dh.CountSize = 10;
dh.recentdatapoints = zeros(1,dh.CountSize);
dh.maxV = zeros(1,dh.NumberOfChannel);
dh.maxI = zeros(1,dh.NumberOfChannel);
dh.MaxAV = 220;
dh.MinAV = 100;
dh.Vdelay = 10;
dh.diff = zeros(1,dh.CountSize);
hisindlearn = zeros(1,length(aindlearn));
n=1;
%use vent ind as groundtruth due to its higher detection accuracy
for k = 1:length(vindlearn(:,1))-1
    
    dh.Vind = vindlearn(k,1)+dh.Vdelay;
    dh.Aind = aindlearn(n);
    if ((dh.Vind - dh.Aind) < dh.MaxAV) && ((dh.Vind - dh.Aind) > dh.MinAV)
        n = n+1;
        atrFound = 't';
    else
        atrFound = 'f';
        dh.Aind =  dh.Vind - dh.MaxAV;
        
    end
    %{
    %see if we find ventricular peak after the atrial peak
    if (vindlearn(k,1) - dh.Aind>dh.MinAV) && (vindlearn(k,1)-dh.Aind)<dh.MaxAV
        ventFound = 't';
    else
        ventFound = 'f';
    end
    %}
    NumToCount = dh.Vind-dh.Aind;

    MaxRepeat = floor(NumToCount/dh.CountSize); 
    hisdata_AV = hisdata(dh.Aind:dh.Vind,:);
    %for each beat,loop to send data and a signal of successful/not atrial detection
    for h = 1:dh.NumberOfChannel
            %start searching each channel by picking largest signal in a chunk of datapoints
        repeat = 0; %how many chunks of datapoints searched
        count = 1;
       % while (atrialFound == 't')  %if found the atrial beat, start his detection
        %start recording 20 datapoints ahead
        while(repeat< MaxRepeat-1)
            while(count < dh.CountSize)
                iHis = hisdata_AV(repeat*dh.CountSize+count,h); %incoming data sample
                dh.recentdatapoints(count) = iHis; %store incoming 10 datapoints
                if count > 1
                    dh.diff(count) = iHis - dh.recentdatapoints(count-1);
                end
                count = count + 1;
            end
            tempmaxV = sum(abs(dh.diff))*sum(abs(dh.recentdatapoints));
            if (tempmaxV > dh.maxV(h))

                dh.maxV(h) = tempmaxV;
                
                [~,I] = max(dh.recentdatapoints);
                dh.maxI(h) = dh.Aind+repeat*dh.CountSize + I;
            end
            repeat = repeat + 1;
            count = 1;
        end
        CountLeft = NumToCount - MaxRepeat*dh.CountSize;
        if (repeat == MaxRepeat-1)
            while(count < CountLeft)
                 iHis = hisdata_AV(MaxRepeat*dh.CountSize+count,h); %incoming data sample
                 dh.recentdatapoints(count) = iHis; %store incoming 10 datapoints
                 if count > 1
                    dh.diff(count) = iHis - dh.recentdatapoints(count-1);
                 end
                count = count + 1;
            end
            tempmaxV = sum(abs(dh.diff))*sum(abs(dh.recentdatapoints));
            if tempmaxV > dh.maxV(h)
                dh.maxV(h) = tempmaxV;
                [~,I] = max(dh.recentdatapoints);
                dh.maxI(h) = dh.Aind+repeat*dh.CountSize + I;
            end
        end

    end
    [~,hchosen] = max(dh.maxV);
    hisI = dh.maxI(hchosen);
    hisV = hisdata(hisI,hchosen);
    HisInfo(k,:) = [hisV,hisI];
    dh.maxV = zeros(1,dh.NumberOfChannel);
    dh.maxI = zeros(1,dh.NumberOfChannel);

end

figure();
hold on;

plot(hisdata(:,1));
plot(hisdata(:,2));
plot(hisdata(:,3));
plot(hisdata(:,4));
plot(hisdata(:,5));
legend('1','2','3','4','5')
plot(HisInfo(:,2),ones(1,length(HisInfo(:,2))),'b*');
