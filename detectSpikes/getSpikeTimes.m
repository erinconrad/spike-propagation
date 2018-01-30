function [gdf,electrodeData] = getSpikeTimes(desiredTimes,dataName,electrodeFile,ptInfo,pwfile,dummyRun)

%{
This is my primary function to detect spikes and output them to a gdf 
with spike times and locations. gdf is an nx2 array, where n is the number
of spikes, the first column has the channel location of the spike and the
2nd column has the time (in seconds) of the spike

You can either call it without arguments, in which case it will use the
hardcoded EEG file listed below to pull from, or you can pass it arguments
from another script.

dummyRun is set to 1 if I am only doing this to generate the electrode
location file, which I do once per patient

It has the ability to use one of two spike detection algorithms: 
- an algorithm by Janca et al. 2014, which detects transient changes in a
signal envelope. I modified this slightly because this begins by
downsampling the data to 200 Hz, which I don't want (higher resolution
important for picking out the order of spikes in the sequence). I also
changed the default filter type from Chebyshev II to FIR because it
appeared to work better with the higher frequency.
- an algorithm that was used in Eric Marsh's lab, which I believe was
written by Camilo Bermudez 7/31/13, and appears to work by finding peaks in
high frequency data and then requires that the peaks have a minimum
amplitude and appropriate duration. I have not modifed this at all yet

Janca et al paper:
Janca, Radek, et al. "Detection of interictal epileptiform discharges using signal envelope distribution modelling: application to epileptic and non-epileptic intracranial recordings." Brain topography 28.1 (2015): 172-183.

%}

%% Parameters to change each time
whichDetector = 1; %1 = modified Janca detector, 2 = Bermudez detector

% If you didn't pass it arguments, then use the ones written here
if nargin == 0
    dataName = 'HUP78_phaseII-Annotations'; % the report you want to look at
    electrodeFile = 'HUP078_T1_19971218_electrode_labels.csv';
    dummyRun = 0;
    
    % These times are in seconds and correspond to the ieeg.org times
    desiredTimes = [267000,268000]; %[267375,267415]; is a seizure % the first and last seconds you want to look at 
    
    
    % the file to write to (will only write if you didn't pass it arguments)
    outFile = 'gdf267000-268000.mat';
    
    % call path file (this defines the locations of various files)
    spikePaths
    fileLocations
    
    % Load pt file (contains seizure times and which electrodes to ignore - like EKG leads)
    ptInfo = loadjson(jsonfile);
    
    
    addpath('/Users/erinconrad/Desktop/residency stuff/R25/actual work/scripts/my scripts/makeChannelStructs/');

end



% Bermudez algorithm parameters
tmul=13; % threshold multiplier
absthresh=300;

%% Parameters that probably don't need to change each time
ignore = 1; % should we ignore any electrodes? This breaks if I say no.



%% Load EEG data info
% calling this with 0 and 0 means I will just get basic info like sampling
% rate and channel labels
data = getiEEGData(dataName,0,0,pwfile);  




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


%% Prep what data I want to look at



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

if dummyRun == 1
    gdf = 0;
    %% make the list of channel locations
    
    electrodeData = chanLocUseGdf(unignoredChLabels,electrodeFile);

elseif dummyRun == 0
    
    % Get the indices I want to look at
    startAndEndIndices = desiredTimes*data.fs;
    indices = startAndEndIndices(1):startAndEndIndices(2);

    % get the data from those indices and channels (ignoring ignored channels)
    data = getiEEGData(dataName,channels,indices,pwfile);

    %% Do cleaning?
    % What should I do?

    %% Run spike detector
    if whichDetector == 1

        % This is the spike detector from Janca et al 2014, edited by me as
        % above
        [out,~,~,~,~,~] = spike_detector_Erin(data.values,data.fs);

        % reorder spikes by time
        [timeSort,I] = sort(out.pos);
        chanSort = out.chan(I);

        % make gdf
        if isempty(out.pos) == 1
            fprintf('No spikes detected\n');
        else
            %fprintf('Detected %d spikes\n',length(out.pos));
            gdf = [chanSort,timeSort];
        end

    elseif whichDetector == 2
         addpath('./SamCode');
        % This calls the Bermudez detector
        % 
        % I have not edited this at all at this point.
        gdf = fspk2(data.values,tmul,absthresh,length(channels),data.fs);


        if isempty(gdf) == 1
            fprintf('No spikes detected\n');
        else
            %fprintf('Detected %d spikes\n',size(gdf,1));
             % put it in seconds
             gdf(:,2) = gdf(:,2)/data.fs;

            out.pos = gdf(:,2); out.chan = gdf(:,1);
        end

    end



    if nargin == 0
        %% Save spike times
        % note that this only has data from the unignored channels
        save([gdfFolder outFile], 'gdf','unignoredChLabels');

        %% Sample plot
        if 1 == 0
        figure
        indicesToPlot = 10000:14000;
        chsToPlot = [68];
        plotSpikeTimes(data,out,indicesToPlot,chsToPlot)
        end
    end
    
    electrodeData = 0;
    
end


end