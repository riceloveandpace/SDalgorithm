function d = singlePeakFinder(i,d)
%This function is called every sample to determine whether the current
%sample is a beat (ventricular/atrial depend on input channel).  To do
%this, the function considers the detection structure, d, and determines
%whether a beat has occured, based on the parameters of d.

%Inputs:
%i - The current sample of the signal.
%d - The detection structure that contains the parameters necessary for
%detection.

%Outputs:
%d.recentBools = a binary vector indicating whether the recent data points
%have been above the threhsold.
%d.beatDelay = delay since the last beat is seen.
%d.PeakInd = time (in samples) of the ventricular peaks.
%d.last_sample_is_beat = Flag indicating whether the last sample was a
%ventricular beat.

%At beginning of function call, assume that no peak has been found.
datapoint = d.recentdatapoints(end);
d.recentBools = [d.recentBools(2:end) (d.flip*datapoint>d.thresh)];

%Only start looking for beats when the current data point
%is greater than the threshold and it has been 30 samples since the
%last ventricle or atrial beat.

%We detect a ventricular beat if the following conditions are met:
%1. If atleast 2/3*v_length of the recent samples have exceeded the
%ventricular threshold.
%2. The delay since the last atrial beat detected is greater than the
%delay since the last atrial falling beat.
%3. It has been at least VV since the last vent beat.
if(sum(d.recentBools)>d.length*2/3 && d.beatDelay >= d.beatFallDelay && d.beatFallDelay > d.VV)
    if (~d.last_sample_is_sig)
        d.beatDelay = 0;
      %  d.beatWeighted = false;
        d.PeakInd = [d.PeakInd, i];
        d.last_sample_is_sig = true;
    end
else
    %If the last sample is V, set the falling edge delay to the sample after.
    if(d.last_sample_is_sig)
        d.beatFallDelay = 0;
        d.last_sample_is_sig = false;
    end
end

end
