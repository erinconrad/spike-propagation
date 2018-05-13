function gdf = portEKG(desiredTimes,dataName,pwfile,tmul,absthresh)

%{
This is my primary function to detect spikes and output them to a gdf 
with spike times and locations. gdf is an nx2 array, where n is the number
of spikes, the first column has the channel location of the spike and the
2nd column has the time (in seconds) of the spike
%}

%% Parameters
whichDetector = 2; %1 = modified Janca detector, 2 = Bermudez detector, 3 = orig Janca


%% Load EEG data info
% calling this with 0 and 0 means I will just get basic info like sampling
% rate and channel labels
data = getiEEGData(dataName,0,0,pwfile);  


%% Find ekg electrode numbers
ekg_electrode_nums = [];
ekg_electrode_labels = {};


 % loop through all channel labels
 for i = 1:length(data.chLabels)   

    %% parsing of channel names (labeled odd in the iEEG)

    % get the name
    origStr = data.chLabels{i};

    if contains(origStr,'EKG') == 1
        ekg_electrode_nums = [ekg_electrode_nums,i];
        ekg_electrode_labels = [ekg_electrode_labels,origStr];
    end

 end

 channels = ekg_electrode_nums;
 fprintf('%d EKG channels found for %s\n',length(ekg_electrode_nums),dataName);


%% Prep what data I want to look at

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

end