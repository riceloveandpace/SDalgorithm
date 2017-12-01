function [a_indices] = atrial_peak_finder(detection,data)
t_blank = 100; l_blank = 100; %number of samples to zero out on each side of a peak

    %find atrial peaks
    a_bool = detection.flip*data>detection.thresh;
    [~,a_indices_r, a_indices_f] = CountPeaks(a_bool,detection.alength);
    a_indices = a_indices_r;
end
