# Signal Detection MATLAB Code Overview

========================

asdf 

### ventDecision.m 

+ Real time decision of peaks on data.


### singlePeakFinder.m 
+ Input: Current Sample, Data Structure (vector of recent samples, delay since last beat, peak index)
+ Output: Data Structure (vector of recent samples, delay, peak index). Updates peak index.
 
### av_detect.m 
+ Event decision, but not real time. Atrial and ventricular algorithm. Looks at the data as a whole.

### atrial_peak_finder.m 
- nonrealtime version for the singlePeakFinder.m

### atrialParamLearning.m 
- parameter learning for atrial and ventricle. Threshold and length of signal. No real time version. 
- Input: Data
- Output: height, flip, a_length
- Subfunctions: 
    - LearnHeightRange(data, minlength)
    - LearnLengthsnew(data)
- Learns parameters for the atrial and ventrical signals.

### His_detect.m 
- N/A

### His_Detection.m 
- N/A

### isDetect_RealTime.m~ 
- Real time code. 
