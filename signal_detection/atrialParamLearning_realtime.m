
function [d] = atrialParamLearning(data,data2,data3,d)

    %learn the appropriate lengths of atria, as well as
    %which should be detected first, for + and - data
    
    
    %learn max value for 3 seconds(for simplicity in matlab, used max function)
    d.max_val = max(data);
    
    %learn height range for 5 seconds
    max_th = d.max_val/2;
    min_th = max_th/10;
    ths = linspace(max_th,min_th,120);
    count = zeros(1,length(ths));
    for i = 1:length(data2)
        s = data2(i);
        for j = 1:length(ths)
            if s > ths(j)
                count(j) = count(j) + 1;
            end
        end
        
        
    end
    temp = diff(diff(count));
    idx = 0;
    for ii = 1:length(temp)
        if temp(ii) < 100
            idx = [idx ii];
        end
          
    end

    d.cutoffs = [ths(idx(end)),ths(idx(2))];
    %[~,af] = max(cutoffs(:,2)-cutoffs(:,1));
    %pick exact thresholds from the middle of the ranges according to mid_finders
    d.thresh = sum(d.cutoffs.*[0.5, 0.5])
    
    d.len = 1;
    
    %learn energy threshold value for 10 seconds
    last_vals = zeros(1,20);
    noiselvls = [];
    count = 1;
    for j = 1:length(data3)
        k = data3(j);
        
        if count < 21
            last_vals(count) = k;
            
        elseif (count > 20) && (k < d.thresh) && (j < d.lastPI + d.VV) 
            noiselvls = [noiselvls mean(abs(last_vals))];
            count = 0;
            last_vals = zeros(1,20); 
        end
        %new beat detected
        if (k > d.thresh) && (j > d.lastPI + d.VV)
            d.lastPI = j;
        end
        count = count + 1;
        
        

        
    end
    
    d.noiselvl = mean(noiselvls);
    %}
    
    
    

    
end
