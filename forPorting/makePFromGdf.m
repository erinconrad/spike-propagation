%% File names
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
%electrodeFile = [electrodeFolder,csvFile];
p1 = genpath(scriptFolder);
addpath(p1);
timeFile = 'desiredTimes.mat';
gdfFolder = 'gdf/';
chLocationsFolder = 'chLocations/';
newptfile = 'ptWithfs.mat';

%% Load file with filenames and run times
load([resultsFolder,timeFile]);

%% Loop through patients, szs, run times
for i = 1:length(pt)
    
    pt(i).fs = 500;
    
end

save([resultsFolder,newptfile],'pt');