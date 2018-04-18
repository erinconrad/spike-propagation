%% make_gdf

% This script takes a patient structure with information about what times
% to look for spikes over, and then this actually detects those spikes and
% outputs them to a gdf file


clear

%% Parameters

% Should I re-run the spike detection and overwrite gdf file if it already
% exists?
overwrite = 0; 

% Should we try to merge the patient structure with an existing, incomplete
% patient structure?
merge = 1;

% Berumudez spike detector parameters
tmul = 13; % the threshold value for the Bermudez spike detector (default is 13)
absthresh = 300; % absolute threshold value for the Bermudez spike detector (default is 300)

%% File names
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
%electrodeFile = [electrodeFolder,csvFile];
p1 = genpath(scriptFolder);
addpath(p1);
timeFile = 'desiredTimes.mat'; 
gdfFolder = [resultsFolder,'gdf/'];
chLocationsFolder = 'chLocations/';
newptfile = 'ptWithfs.mat';

%% Load file with filenames and run times
if merge == 1 && exist([resultsFolder,'ptStructs/',newptfile],'file') ~= 0
    load([resultsFolder,'ptStructs/',newptfile]);
else
    load([resultsFolder,'ptStructs/',timeFile]);
end

%% Loop through patients, szs, run times
for i = 1:length(pt)
    
    dataName =  pt(i).ieeg_name;
    if isempty(dataName) == 1
        continue
    end
    
    electrodeFile = pt(i).electrode_labels;
    if isempty(electrodeFile) ==1
        continue
    end
    
    ignoreElectrodes = pt(i).ignore_electrodes;
    
    for j = 1:length(pt(i).sz)
        if isfield(pt(i).sz(j),'runTimes') == 0
            continue
        end
        for k = 1:size(pt(i).sz(j).runTimes,1)
            
            % Add a button push to the desmond file (for the purpose of
            % restarting the program if it crashes due to random server
            % error)
            buttonpush = datestr(now,'yyyy-mm-dd HH:MM:SS');
            allwrite = [buttonpush,'\n',sprintf('Patient %s seizure %d chunk %d\n',...
                dataName,j,k)];
            fid = fopen('/tmp/desmond.txt','wt');
            fprintf(fid,allwrite);
            fclose(fid);
            
            if exist([gdfFolder,pt(i).name,'/',pt(i).sz(j).chunkFiles{k}],'file') ~= 0
                if overwrite == 0
                    fprintf('File %s already found, skipping\n',[pt(i).sz(j).chunkFiles{k}]);
                    continue
                end
            end
            
            desiredTimes = [pt(i).sz(j).runTimes(k,:)];
            [gdf,electrodeData,fs] = portGetSpikes(desiredTimes,dataName,...
                [electrodeFolder,electrodeFile],ignoreElectrodes,pwfile,tmul,absthresh);
            
            % Add fs and electode data to the patient structure
            pt(i).fs = fs;
            pt(i).electrodeData = electrodeData;
            
            % Save gdf file
            save([gdfFolder,pt(i).name,'/',pt(i).sz(j).chunkFiles{k}],'gdf');
            
            % Save chLocations file
            save([electrodeFolder,pt(i).chLocationFile],'electrodeData');
            
           
            % Resave pt file now that I have fs
            save([resultsFolder,'ptStructs/',newptfile],'pt');
            
        end
    end
end



% Make a new document if I make it here
fid2 = fopen('/tmp/ok.txt','wt');
fprintf(fid2,'Done\n');
%fflush(fid2);
fclose(fid2);
