clear

%% Parameters

% Should I re-run the spike detection and overwrite gdf file if it already
% exists?
overwrite = 0; 

%% File names
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
%electrodeFile = [electrodeFolder,csvFile];
p1 = genpath(scriptFolder);
addpath(p1);
timeFile = 'finalPt.mat';
gdfFolder = 'gdf/';
chLocationsFolder = 'chLocations/';
newptfile = 'ptWithfs.mat';

%% Load file with filenames and run times
load([resultsFolder,timeFile]);

%% Loop through patients, szs, run times
for i = 1:length(pt)
    
    
    dataName = pt(i).ieeg_name;
    if isempty(dataName) == 1
        continue
    end
    electrodeFile = pt(i).electrode_labels;
    if strcmp(electrodeFile,'??') ==1
        continue
    end
    
    ignoreElectrodes = pt(i).ignore_electrodes;
    
    data = getiEEGData(dataName,0,0,pwfile);  
    fs = data.fs;

    if fs ~= 512
        fprintf('Warning, patient %d (%s) has an fs of %d\n',i,pt(i).name,fs);
    else
        fprintf('It is ok, patient %d (%s) had correct fs\n',i,pt(i).name);
    end
        
    
end

