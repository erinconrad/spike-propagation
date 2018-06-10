%% getElectrodeData
clear

overwrite =  1;

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
    
    ignoreElectrodes = pt(i).ignore_electrodes;
    
    for j = 1:length(pt(i).sz)
        if isfield(pt(i).sz(j),'runTimes') == 0
            continue
        end
        
        %% Load EEG data info
        % calling this with 0 and 0 means I will just get basic info like sampling
        % rate and channel labels
        data = getiEEGData(dataName,0,0,pwfile);  
        
        %% Ignore certain electrodes (EKG, etc.)
        % Make cell of channel names
        chNames = cell(length(data.chLabels),1);

        % initialize which channels to ignore
        chIgnore = zeros(length(data.chLabels),1);
        nchan = length(data.chLabels);

         % loop through all channel labels
         for k = 1:length(data.chLabels)   

            %% parsing of channel names (labeled odd in the iEEG)

            % get the name
            origStr = data.chLabels{k};

            % split it with spaces
            C = strsplit(origStr,' ');

            % I would expect all of the names to start with EEG
            if strcmp(C{1},'EEG') == 0
                fprintf('Warning, there is something weird in the channel labels for channel %d in patient %s\n',i,dataName);
                C = C{1};

            else
                C = strrep(origStr,[C{1},' '],'');

                % Remove -Ref
                D = strsplit(C,'-');

                C = strrep(C,['-',D{2}],'');
            end


            % Remove space if present
            C = strrep(C,' ','');

            % Get the numbers
            numIdx = regexp(C,'\d');

            if isempty(numIdx) == 0
                if strcmp(C(numIdx(1)),'0') == 1
                    C(numIdx(1)) = [];
                end
            end

            % Final channel name
            chName = C;
            chNames{k} = chName;

         end
         
         for x = 1:length(data.chLabels)
            chName = chNames{x};
            for y = 1:length(ignoreElectrodes)
                if strcmp(chName,ignoreElectrodes{y}) == 1
                    chIgnore(x) = 1;
                end
            end
    
         end
         
         channels = find(chIgnore == 0);
         unignoredChLabels = chNames(channels);
         nchan = length(channels);
         electrodeData = chanLocUseGdf(unignoredChLabels,[electrodeFolder,'/',electrodeFile]);


        % Add fs and electode data to the patient structure
        pt(i).fs = data.fs;
        pt(i).electrodeData = electrodeData;
        pt(i).channels = channels;


        % Save chLocations file
        save([electrodeFolder,pt(i).chLocationFile],'electrodeData');

            
    end
            
end

save([resultsFolder,'ptStructs/',newptfile],'pt');

