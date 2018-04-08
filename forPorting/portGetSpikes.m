function [gdf,electrodeData,fs] = portGetSpikes(desiredTimes,dataName,...
    electrodeFile,ignoreElectrodes,pwfile)

%{
This is my primary function to detect spikes and output them to a gdf 
with spike times and locations. gdf is an nx2 array, where n is the number
of spikes, the first column has the channel location of the spike and the
2nd column has the time (in seconds) of the spike
%}

%% Parameters
whichDetector = 2; %1 = modified Janca detector, 2 = Bermudez detector, 3 = orig Janca

setChLimits = 0;
multiChLimit = 0.8; % I will throw out spikes that occur in >80% of channels at the same time
multiChTime = .025; % The time period over which spikes need to occur across multiple channels to toss

% Bermudez algorithm parameters
tmul=13; % threshold multiplier
absthresh=300;


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

 end

   
for i = 1:length(data.chLabels)
    chName = chNames{i};
   
    for j = 1:length(ignoreElectrodes)
        if strcmp(chName,ignoreElectrodes{j}) == 1
            chIgnore(i) = 1;
        end
    end
    
end

channels = find(chIgnore == 0);
unignoredChLabels = chNames(channels);
nchan = length(channels);



%% Prep what data I want to look at

% get the channels I want to look at. Also make unignoredChLabels, which is
% the length of the new channel array and keeps track of the channel
% identities of these newly indexed channels



%% make the list of channel locations
electrodeData = chanLocUseGdf(unignoredChLabels,electrodeFile);

% Get the indices I want to look at
startAndEndIndices = desiredTimes*data.fs;
indices = startAndEndIndices(1):startAndEndIndices(2);

% get the data from those indices and channels (ignoring ignored channels)
data = getiEEGData(dataName,channels,indices,pwfile);


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
    
fs = data.fs;

end