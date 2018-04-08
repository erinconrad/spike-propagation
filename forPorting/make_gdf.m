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
timeFile = 'desiredTimes.mat';
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
    
    for j = 1:length(pt(i).sz)
        if isfield(pt(i).sz(j),'runTimes') == 0
            continue
        end
        for k = 1:size(pt(i).sz(j).runTimes,1)
            
            if exist([pt(i).sz(j).chunkFiles{k}],'file') ~= 0
                if overwrite == 0
                    continue
                end
            end
            
            desiredTimes = [pt(i).sz(j).runTimes(k,:)];
            [gdf,electrodeData,fs] = portGetSpikes(desiredTimes,dataName,...
                electrodeFile,ignoreElectrodes,pwfile);
            
            pt(i).fs = fs;
            
            % Save gdf file
            save([pt(i).sz(j).chunkFiles{k}],'gdf');
            
            % Save chLocations file
            save([pt(i).chLocationFile],'electrodeData');
            
           
            
            
        end
    end
end

% Resave pt file now that I have fs
save([resultsFolder,newptfile],'pt');
