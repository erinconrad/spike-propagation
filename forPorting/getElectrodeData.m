%% getElectrodeData
clear

% 2 means I ignore channels if and only if they don't have a location in 
% the electrode csv file
useErinIgnore = 2; 

overwrite =  0;

%% File names
newptfile = 'ptWithElectrodeData.mat';
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
%electrodeFile = [electrodeFolder,csvFile];
p1 = genpath(scriptFolder);
addpath(p1);
gdfFolder = [resultsFolder,'gdf/'];
ptInfo = loadjson(jsonfile);
ptnames = fieldnames(ptInfo.PATIENTS);
timeFile = 'desiredTimes.mat';

%% Load file with filenames and run times
if overwrite == 0 && exist([resultsFolder,'ptStructs/',newptfile],'file') ~= 0
    load([resultsFolder,'ptStructs/',newptfile]);
else
    load([resultsFolder,'ptStructs/',timeFile]);
end


%% Loop through patients, szs, run times
for i = 1:length(pt)
    fprintf('Doing patient %s\n',pt(i).name);
    dataName =  pt(i).ieeg_name;
    if isempty(dataName) == 1
        continue
    end
    
    electrodeFile = pt(i).electrode_labels;
    if isempty(electrodeFile) ==1
        continue
    end
    
    if overwrite == 0 && isfield(pt(i),'electrodeData') == 1 && isempty(pt(i).electrodeData) == 0
        fprintf('Got electrode data for patient %s, skipping...\n',pt(i).name);
        continue 
    end
    
    if useErinIgnore == 0
        ignoreElectrodes = pt(i).ignore_electrodes;
    elseif useErinIgnore == 1
        ignoreElectrodes = pt(i).erin_ignore;
    elseif useErinIgnore == 2
        
    end
    
    for j = 1:length(pt(i).sz)
        if isfield(pt(i).sz(j),'runTimes') == 0
            continue
        end
        
        %% Load EEG data info
        % calling this with 0 and 0 means I will just get basic info like sampling
        % rate and channel labels
        data = getiEEGData(dataName,0,0,pwfile);  
        
        
        %% Ignore certain electrodes (EKG, etc.)
        
        chLabels = data.chLabels;
        chLabelsParsed = chLabels;
        nchan = length(chLabelsParsed);
        chIgnore = zeros(nchan,1);
        
        % Get the channels labels from ieeg and then parse them to be
        % easier to read and compare to my names
        for ch = 1:length(chLabels)
           chLabelsParsed{ch} = chParser(chLabels{ch}); 
        end
        
        % Ignore any channels that do not have locations in the electrode
        % location file. These are the ONLY channels I will ignore.
        if useErinIgnore == 2
            ignoreElectrodes = findChsToIgnore(pt,i,chLabelsParsed);
        end
       
                 
        foundIgnoredChs = zeros(length(ignoreElectrodes),1);
        
        % Find channels that are equal to my ignored channels
        for x = 1:length(chLabelsParsed)
            chName = chLabelsParsed{x};
            for y = 1:length(ignoreElectrodes)
                if strcmp(chName,ignoreElectrodes{y}) == 1
                    chIgnore(x) = 1;
                    foundIgnoredChs(y) = 1;
                end
            end
        end
         
        % Give a warning if I have remaining ignoreElectrodes that I could
        % not identify amongst the channels
        if sum(foundIgnoredChs) ~= length(ignoreElectrodes)
            unfoundIdx = find(foundIgnoredChs==0);
            unfound = ignoreElectrodes{unfoundIdx};
            fprintf('Warning, could not find ignored channel %s\n',unfound);
        end
        
        channels = find(chIgnore == 0);
        unignoredChLabels = chLabelsParsed(channels);
        nchan = length(channels);
        electrodeData = chanLocUseGdf(unignoredChLabels,[electrodeFolder,'/',electrodeFile]);
        electrodeData.allLabels = data.chLabels;
        
        % Add fs and electode data to the patient structure
        pt(i).fs = data.fs;
        pt(i).electrodeData = electrodeData;
        pt(i).channels = channels;

     

        % Save chLocations file
        save([electrodeFolder,pt(i).chLocationFile],'electrodeData');

            
    end
            
end

save([resultsFolder,'ptStructs/',newptfile],'pt');

