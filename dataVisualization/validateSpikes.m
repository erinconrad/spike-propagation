function validateSpikes(pt,whichPt,tmul,absthresh,whichTime,whichChs)

%% Parameters
% how much time to detect
duration = 600; %default 600, less than this will screw up sensitivity of detector

% how much time to check for spikes
spike_dur =  1; %1 s before, 1 s ater


%% Initialize
% Get paths
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;

dataName = Patient(pt).ieeg_name;
electrodeFile = Patient(pt).electrode_labels;
ptname = Patient(pt).name;


end