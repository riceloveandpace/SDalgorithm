# Overview of code in this folder

ventDecision.m - real time decision of peaks on data.

singlePeakFinder.m - subfunction. Called each data point. Living on the node. Detect if it is a peak or not.
 Input: Current Sample, Data Structure (vector of recent samples, delay since last beat, peak index)
 Output: Data Structure (vector of recent samples, delay, peak index). Updates peak index.
 
av_detect.m - Event decision, but not real time. Atrial and ventricular algorithm. Looks at the data as a whole.

atrial_peak_finder.m - nonrealtime version for the singlePeakFinder.m

atrialParamLearning.m - parameter learning for atrial and ventricle. Threshold and length of signal. No real time version. 

His_detect.m - N/A

His_Detection.m - N/A

HisDetect_RealTime.m~ - Real time code. 
