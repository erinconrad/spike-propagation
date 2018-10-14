function pt = updateGDFNewTimes(whichPt)

%% File names
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
%electrodeFile = [electrodeFolder,csvFile];
p1 = genpath(scriptFolder);
addpath(p1);
timeFile = 'ptWithElectrodeData.mat'; 
gdfFolder = [resultsFolder,'gdf/'];
newptfile = 'ptPostGDF.mat';

%% Load file with filenames and run times
ptGDF = load([resultsFolder,'ptStructs/',newptfile]);

ptTimes = load([resultsFolder,'ptStructs/',timeFile]);


ptGDF.pt(whichPt) = ptTimes.pt(whichPt);

pt = ptGDF.pt;

save([resultsFolder,'ptStructs/',newptfile],'pt');
