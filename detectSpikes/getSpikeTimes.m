function [gdf,electrodeData,extraOutput] = getSpikeTimes(desiredTimes,json_name,dataName,...
    electrodeFile,ptInfo,pwfile,dummyRun,vanleer,vtime,outputData,keepEKG,ignore,funnyname,tmul,absthresh)

%{
This is my primary function to detect spikes and output them to a gdf 
with spike times and locations. gdf is an nx2 array, where n is the number
of spikes, the first column has the channel location of the spike and the
2nd column has the time (in seconds) of the spike

You can either call it without arguments, in which case it will use the
hardcoded EEG file listed below to pull from, or you can pass it arguments
from another script.

dummyRun is set to 1 if I am only doing this to generate the electrode
location file, which I do once per patient. The reason I do this in the
spike detector script is because I want to ensure that I consistently track
which electrodes I am not running the spike detector on, and so in this
script I look up which electrodes to ignore as listed in the json file, and
then I store the identities of those electrodes in my electrodeData file
and I ignore them when I run the spike detector.

It has the ability to use one of two spike detection algorithms: 
- an algorithm by Janca et al. 2014, which detects transient changes in a
signal envelope.
- an algorithm that was used in Eric Marsh's lab, which I believe was
written by Camilo Bermudez 7/31/13, and appears to work by finding peaks in
high frequency data and then requires that the peaks have a minimum
amplitude and appropriate duration. I have not modifed this at all yet
    
Other parameters:
- vanleer: 1 if I'm doing the vanleer method
- vtime: time for the vanleer method
- outputData: 1 if I want to output the EEG data along with the spikes 
- keepEKG: 1 if I just want to pull out the EKG channel
- ignore: 1 if I ignore the channels listed in the json file, 0 if I ignore
nothing
- funnyname: 0 if it's one of the standard-named HUP patients who I can
look up in the json file. 1 if it's another type of dataset

Janca et al paper:
Janca, Radek, et al. "Detection of interictal epileptiform discharges using signal envelope distribution modelling: application to epileptic and non-epileptic intracranial recordings." Brain topography 28.1 (2015): 172-183.

%}

%% Parameters
whichDetector = 2; %1 = modified Janca detector, 2 = Bermudez detector, 3 = orig Janca
timeToLookForPeak = .1; % look 100 ms before and 100 ms after the detected spike to find the peak 

setChLimits = 0;
multiChLimit = 0.8; % I will throw out spikes that occur in >80% of channels at the same time
multiChTime = .025; % The time period over which spikes need to occur across multiple channels to toss

% Bermudez algorithm parameters
%tmul=10; % threshold multiplier, default is 13
%absthresh=300;

% If you didn't pass it arguments, then use the ones written here
if nargin == 0
    dataName = 'HUP78_phaseII-Annotations'; % the report you want to look at
    electrodeFile = 'HUP078_T1_19971218_electrode_labels.csv';
    dummyRun = 0;
    
    % These times are in seconds and correspond to the ieeg.org times
    desiredTimes = [267375,267415]; %[267375,267415]; is a seizure % the first and last seconds you want to look at 
    
    
    % the file to write to (will only write if you didn't pass it arguments)
    outFile = 'gdf267000-268000.mat';
    
    % call path file (this defines the locations of various files)
    spikePaths
    [electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
    
    % Load pt file (contains seizure times and which electrodes to ignore - like EKG leads)
    ptInfo = loadjson(jsonfile);
    vanleer = 0;
    outputData = 1;
    keepEKG = 0;
    addpath('/Users/erinconrad/Desktop/residency stuff/R25/actual work/scripts/my scripts/makeChannelStructs/');
    ignore = 0; % should we ignore any electrodes? 
end

ptname = json_name;

%% Load EEG data info
% calling this with 0 and 0 means I will just get basic info like sampling
% rate and channel labels
data = getiEEGData(dataName,0,0,pwfile);  


%% Select correct patient in order to get ignore electrode info from json file
% requires parsing because the name of the patient in json file is
% different from iEEG
% only look at the part of the name before the _




%% Get electrodes to ignore
% Make cell of channel names
chNames = cell(length(data.chLabels),1);

% initialize which channels to ignore
chIgnore = zeros(length(data.chLabels),1);
nchan = length(data.chLabels);
channels = 1:nchan;

     % loop through all channel labels
     for i = 1:length(data.chLabels)   
         
         % Don't need to do this for non-HUP datasets as I don't have them
         % in the json file
         if funnyname == 0
            

             %% parsing of channel names (labeled odd in the iEEG)

            % get the name
            origStr = data.chLabels{i};

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
            chNames{i} = chName;
         else
             origStr = data.chLabels{i};
             chName = origStr;
             chNames{i} = chName;
         end
    end

if ignore == 1
% now we can look up the name of the patient in the json file format to see
% which electrodes to ignore
ignoreElectrodes = ptInfo.PATIENTS.(ptname).IGNORE_ELECTRODES;

   
    for i = 1:length(data.chLabels)
        chName = chNames{i};
        %% ignore certain electrodes as suggested in the json file
        if keepEKG == 0
        
        % Ignore certain electrodes
            for j = 1:length(ignoreElectrodes)
                if strcmp(chName,ignoreElectrodes{j}) == 1
                    chIgnore(i) = 1;
                end
            end
        elseif keepEKG == 1
            % only keep the EKG
           if contains(chName,'EKG') == 1
               chIgnore(i) = 0;
           else
               chIgnore(i) = 1;
           end
       
        end

    end
    
    channels = find(chIgnore == 0);
    unignoredChLabels = chNames(channels);
    nchan = length(channels);
end


%% Prep what data I want to look at

% get the channels I want to look at. Also make unignoredChLabels, which is
% the length of the new channel array and keeps track of the channel
% identities of these newly indexed channels


if dummyRun == 1
    
    % Don't need gdf for dummy run
    gdf = 0;
    %% make the list of channel locations
    electrodeData = chanLocUseGdf(unignoredChLabels,electrodeFile);

elseif dummyRun == 0
    
    % Don't need electrode Data for non-dummy run as already calculated it
    electrodeData = 0;
    
    % Get the indices I want to look at
    startAndEndIndices = desiredTimes*data.fs;
    indices = startAndEndIndices(1):startAndEndIndices(2);

    % get the data from those indices and channels (ignoring ignored channels)
    data = getiEEGData(dataName,channels,indices,pwfile);
    
    %% Remove channels that are just nans and zeros
    if ignore == 0
        allsum = zeros(nchan,1);
        for i = 1:nchan
        allsum(i) = nansum(abs(data.values(:,i)));
        end
        nanidx = (allsum == 0);
        data.values = data.values(:,~nanidx);
        channels = find(nanidx == 0);
        unignoredChLabels = chNames(channels);
        nchan = length(channels);
    end

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
            fprintf('Warning: No spikes detected\n');
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
            fprintf('Detected %d spikes\n',size(gdf,1));
             % put it in seconds
             gdf(:,2) = gdf(:,2)/data.fs;

            out.pos = gdf(:,2); out.chan = gdf(:,1);
        end
        
    elseif whichDetector == 3
        %this is the unedited Janca detector
        
        [out,~,~,~,~,~] = spike_detector_hilbert_v16_nodownsample(data.values,data.fs,'-h 60');

        % reorder spikes by time
        [timeSort,I] = sort(out.pos);
        chanSort = out.chan(I);

        % make gdf
        if isempty(out.pos) == 1
            fprintf('Warning: No spikes detected\n');
        else
            fprintf('Detected %d spikes\n',length(out.pos));
            gdf = [chanSort,timeSort];
        end

    end
    
    %% Re-define the spike time as the peak
    %{
    
    gdfnew = gdf;
    % convert times to indices
    gdfIdx = round([gdf(:,1),gdf(:,2)*data.fs]);
    timeToLookForPeakIdx = round(timeToLookForPeak*data.fs);
    
    % Loop through all spikes
    for i = 1:size(gdf,1)
        
       % define the period over which to search for the peak relative to
       % the detected spike time (50 ms total)
       spikePeriod = [gdfIdx(i,2)-timeToLookForPeakIdx:gdfIdx(i,2)+timeToLookForPeakIdx];
       
       % make sure it's not less than 1 or more than the max of the data
       spikePeriod = spikePeriod(spikePeriod>0);
       spikePeriod = spikePeriod(spikePeriod<size(data.values,1));
      
       % get the channel of the spike
       ch = gdf(i,1);
       
       % Get the index of the maximum
       [~,I] = max(data.values(spikePeriod,ch));
       
       % Get the index relative to the whole block
       newSpikeTimeIdx = spikePeriod(1)-1+I;
       
       % put it back into seconds
       newSpikeTime = newSpikeTimeIdx/data.fs;
       gdfnew(i,2) = newSpikeTime;
        
    end
    gdf=gdfnew;
    
    
    %% reorder gdf
    [timesort,I] = sort(gdf(:,2));
    chanSort = gdf(I,1);
    gdf = [chanSort,timesort];
    %}
    
    %% Toss spikes that occur across too high a percentage of channels at the same time
    
    if setChLimits == 1 && isempty(gdf) == 0
        newgdf =  gdf;
        tooManyChs = [];
        i = 1; % start with i = 1 (the first spike)
        % Loop through the spikes (we are going to variably move through the spikes)
        while 1
            scount = 1;

            % if i is the last spike, break
            if i == size(newgdf,1)
                break
            end

           % for each spike, loop through subsequent spikes to see how many
           % there are across multiple channels occuring at the same time
           for j = i+1:size(newgdf,1)

               % Get the time difference between the first spike and the last
               % spike
               tdiff = newgdf(j,2)-newgdf(i,2);

               % If the difference is small enough, increase the count
               if tdiff < multiChTime
                   scount = scount + 1;
               end

               % if the total number of channels spiking in this very close
               % proximity is >80% of the total number of channels
               if scount > multiChLimit*nchan

                    % Remove these spikes
                    tooManyChs = [tooManyChs;newgdf(i:j,:)];
                    newgdf(i:j,:) = [];

                    % Break the inner for loop and move to the next spike
                    i = j;
                    break
               end

           end

           % advance the index of the spike
           i = i+1;
        end
        fprintf('%d total spikes detected\n',size(gdf,1));
        fprintf('Percentage of spikes discarded for being across too many channels: %1.1f\n',...
            size(tooManyChs,1)/size(gdf,1)*100)
        fprintf('%d spikes remain\n',size(newgdf,1));
        gdf = newgdf;
    end
    
    %% For vanleer
    if vanleer == 1
        % If doing the vanleer approach, input the gdf into a separate
        % function to get delay and rms info
        gdf = vMakeSegments(gdf,data.values,data.fs,vtime);
        
    end

    if nargin == 0
        %% Save spike times
        % note that this only has data from the unignored channels
        save([gdfFolder outFile], 'gdf','unignoredChLabels');
    end
    %% Sample plot
    if 1 == 0
    figure
    indicesToPlot = 10000:14000;
    chsToPlot = [10,30,50,68,80];
    plotSpikeTimes(data,out,indicesToPlot,chsToPlot)
    end
    
    
   
    
end

if outputData == 1
   extraOutput = {data.values,unignoredChLabels,tmul};
else
   extraOutput = 0; 
end

end