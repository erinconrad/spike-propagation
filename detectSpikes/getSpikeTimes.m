clear

%{
This is my primary function to detect spikes and output them to a gdf file
with spike times and locations. It uses the spike detection algorithm by
Janca et al. 2014, which detects transient changes in a signal envelope.

Janca et al paper:
Janca, Radek, et al. "Detection of interictal epileptiform discharges using signal envelope distribution modelling: application to epileptic and non-epileptic intracranial recordings." Brain topography 28.1 (2015): 172-183.

%}

%% Parameters to change each time
dataName = 'HUP78_phaseII-Annotations'; % the report you want to look at
desiredTimes = [267000,268000]; %[267329,267419]; % the first and last seconds you want to look at 

%% Parameters that probably don't need to change each time
ignore = 1; % should we ignore any electrodes? This breaks if I say no.

% call path thing
spikePaths

% file with names of channels to ignore
%% Load EEG data info

% calling this with 0 and 0 means I will just get basic info like sampling
% rate and channel labels
data = getiEEGData(dataName,0,0);  


%% Load pt file
ptInfo = loadjson(jsonfile);

%% Select correct patient
% requires parsing because the name of the patient in json file is
% different from iEEG

% only look at the part of the name before the _
C = strsplit(dataName,'_');
ptname = C{1};

% Find the indices in the name with the numbers
ptnum = regexp(ptname,'\d');

% if the number is only 2 digits, pad it with a leading zero
if length(ptnum) == 2
    ptname = [ptname(1:ptnum(1)-1),'0',ptname(ptnum(1):end)];
elseif length(ptnum) < 2
    fprintf('The name of the patient is unexpected\n');
end

% now we can look up the name of the patient in the json file format to see
% which electrodes to ignore
ignoreElectrodes = ptInfo.PATIENTS.(ptname).IGNORE_ELECTRODES;


%% Get electrodes to ignore
% Make cell of channel names
chNames = cell(length(data.chLabels),1);

% initialize which channels to ignore
chIgnore = zeros(length(data.chLabels),1);

if ignore == 1

    % loop through all channel labels
    for i = 1:length(data.chLabels)
        %% parsing of channel names (labeled odd in the iEEG)
        
        % get the name
        origStr = data.chLabels{i};

        % split it with spaces
        C = strsplit(origStr,' ');

        % I would expect all of the names to start with EEG
        if strcmp(C{1},'EEG') == 0
            fprintf('Warning, there is something weird in the channel labels for channel %d\n',i);
        end

        % Remove leading zero from the number
        if strcmp(C{3}(1),'0') == 1
            endOfChan = C{3}(2:end);
        else
            endOfChan = C{3};
        end

        % Remove -Ref
        D = strsplit(endOfChan,'-');

        % Final channel name
        chName = [C{2},D{1}];
        chNames{i} = chName;

        %% ignore certain electrodes as suggested in the json file
        
        % Ignore certain electrodes
        for j = 1:length(ignoreElectrodes)
            if strcmp(chName,ignoreElectrodes{j}) == 1
                chIgnore(i) = 1;
            end
        end

    end

end

if 1== 1

%% Prep what data I want to look at

% Get the indices I want to look at
startAndEndIndices = desiredTimes*data.fs;
indices = startAndEndIndices(1):startAndEndIndices(2);

% get the channels I want to look at. Also make unignoredChLabels, which is
% the length of the new channel array and keeps track of the channel
% identities of these newly indexed channels
channels = find(chIgnore == 0);
unignoredChLabels = cell(length(channels),1);
j = 0;
for i = 1:length(data.chLabels)
    if chIgnore(i) == 0
        j = j+1;
        unignoredChLabels{j} = chNames{i};
    end
end

% get the data from those indices and channels (ignoring ignored channels)
data = getiEEGData(dataName,channels,indices);

%% Do cleaning
% What should I do?

%% Run spike detector

% This is the spike detector from Janca et al 2014
[out,MARKER,envelope,background,discharges,envelope_pdf] = ...
    spike_detector_hilbert_v16_byISARG(data.values,data.fs);

%% reorder spikes by time
[timeSort,I] = sort(out.pos);
chanSort = out.chan(I);


%% Save spike times
% note that this only has data from the unignored channels
fn = data.name;
if isempty(out.pos) == 1
   fprintf('No spikes detected\n');
else
   fprintf('Detected %d spikes\n',length(out.pos));
   gdf = [chanSort,timeSort];
end
save([gdfFolder ptname '_gdf.mat'], 'gdf','unignoredChLabels');

%% Sample plot
if 1 == 0
indicesToPlot = 6000:10000;
chsToPlot = [10,60];
plotSpikeTimes(data,out,indicesToPlot,chsToPlot)
end

end