% currently: 200 datapoints overlap between each labeling frames; can
% increase if necessary; the overlapped portion of the signal has a
% different color to avoid labelling the same signal twice

% Keyboard input rules:
% red for ventricle, green for atrial, blue for his, magneto for repo;
% n is for next (next string of signal); end is for stop everything;

% after stopped three files will be saved into .mat
% the saved .mat file is an nx2 matrix, each row signifies start - end of a
% peak
function labeler
dataFolder = '/Users/yujunchen/文稿/rice/seniordesign/Data/EPStudys/'; % input
filename = 'ep1TerminatingTachy.mat';%'ep1SVT.mat';%'ep1SR.mat'; % input
initpt = 1; % starting beat; input
in1 = input('Enter channel: "his","atr", or "vent"','s');
if in1 == 'his'
    chnl = [10,11,12,13,14,15];
elseif in1 == 'atr'
    chnl = 16;
elseif in1 == 'vent'
    chnl = 17;
end
overlap = 200; % specify overlay; input
labelergui(dataFolder,filename,initpt,chnl,overlap)
end

function labelergui(dataFolder,filename,initpt,chnl,overlap)
temp = load([dataFolder,filename]);
dat = cell2mat(struct2cell(temp)); % read data into a matrix
if length(chnl) == 1
    sgnl = dat(initpt:end,chnl); % go to the current channel and read out the input
else
    sgnl = dat(initpt:end,chnl(1):chnl(end));
end
cont = 1; idx = 0;
ven = []; atr = []; his = []; rep = [];
% the vector of time stamp to be saved; can be expanded in the function manually

figure(); hold on
chunk = sgnl(1:1000,:);
for i = 1:size(chunk,2)
    plot(chunk(:,i))
end
H = uicontrol();
title(strcat('filename ',', channel:',num2str(chnl),', starting at:',num2str(initpt)))
while (cont == 1)
    % click rule: MUST receive a pair of inputs as start and end point of
    % signal; if one is not included, click on the outside;
    % each different peak is associated and labeled with diff color
    clck = input('Enter color: (r - vent, b - his, m - repo, or g - atr), or hit n to go to the next frame, or hit end to stop the program','s');
    if clck == 'n'
        idx = idx + 1;
        chunk = sgnl(idx*overlap+1:idx*overlap+1000);
        close;
        figure(); plot(chunk,'b'); hold on; plot(chunk(1:(1000-overlap)),'r'); H = uicontrol();
        if length(chnl) == 5
            title(strcat('filename ',', channel: his, starting at:',num2str(initpt+idx*overlap)))
        elseif chnl == 15
            title(strcat('filename ',', channel: atrial, starting at:',num2str(initpt+idx*overlap)))
        elseif chnl == 16
            title(strcat('filename ',', channel: ventricle, starting at:',num2str(initpt+idx*overlap)))            
        end
    elseif clck == 'end'
        cont = 0;
    else 
        allx = [0 0];
        for q = 1:2
            [x,y] = ginput(1);
            x = round(x); y = chunk(x);
            pt = plot(x,y,'*','color',clck);
            satisfy = input('hit enter if you like it and hit z if you do not','s');
            if isempty(satisfy) % then save in this step
                x = initpt-1 + x + idx*(1000-overlap); % add the offset back to get real timestamp
                allx = allx + x*[1,0]*(q==1)+x*[0,1]*(q==2);
            elseif satisfy == 'z'
                disp('Redraw: this point is not saved');
                set(pt,'visible','off')
                ridx = 1;
                while (ridx == 1)
                    [x,y] = ginput(1);
                    x = round(x); y = chunk(x);
                    pt = plot(x,y,'*','color',clck);
                    satisfy = input('hit enter if you like it and hit z if you do not','s');
                    if isempty(satisfy)
                        x = initpt-1 + x + idx*(1000-overlap);
                        allx = allx + x*[1,0]*(q==1)+x*[0,1]*(q==2);
                        ridx = 0;
                    elseif satisfy == 'z'
                        disp('Redraw: this point is not save');
                        set(pt,'visible','off')
                    end
                end
            end
        end
        if clck == 'r'; ven = [ven;allx];
        elseif clck == 'g'; atr = [atr;allx];
        elseif clck == 'b'; his = [his;allx];
        elseif clck == 'm'; rep = [rep;allx];
        end
        
    end
end
save(strcat(filename,num2str(initpt),'_ventricular_timestamp.mat'),'ven')
save(strcat(filename,num2str(initpt),'_his_timestamp.mat'),'his')
save(strcat(filename,num2str(initpt),'_atrial_timestamp.mat'),'atr')
save(strcat(filename,num2str(initpt),'_repolarization_timestamp.mat'),'rep')
end
