whichPt = 31;
dataName = 'Study 029';

%% Get important file paths
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile,other] = fileLocations;
p1 = genpath(scriptFolder);
addpath(p1);
addpath(genpath(other.ieeg));

% Need to check that channels right
times = pt(whichPt).runTimes;
channels = pt(whichPt).channels;
startAndEndIndices = times*pt(whichPt).fs;
indices = startAndEndIndices(1):startAndEndIndices(2);

data = getiEEGData(dataName,channels,indices,pwfile);
pt(whichPt).eeg_data = data;