function gdf = portEKG(desiredTimes,dataName,pwfile,tmul,absthresh,whichDetector,fs,allLabels)

%{
This is my primary function to detect spikes and output them to a gdf 
with spike times and locations. gdf is an nx2 array, where n is the number
of spikes, the first column has the channel location of the spike and the
2nd column has the time (in seconds) of the spike
%}

%% Parameters


%% Find ekg electrode numbers
ekg_electrode_nums = [];
ekg_electrode_labels = {};


 % loop through all channel labels
 for i = 1:length(allLabels)   

    %% parsing of channel names (labeled odd in the iEEG)

    % get the name
    origStr = allLabels{i};

    if contains(origStr,'EKG') == 1
        ekg_electrode_nums = [ekg_electrode_nums,i];
        ekg_electrode_labels = [ekg_electrode_labels,origStr];
    end

 end

 channels = ekg_electrode_nums;
 fprintf('%d EKG channels found for %s\n',length(ekg_electrode_nums),dataName);


%% Prep what data I want to look at

% Get the indices I want to look at
startAndEndIndices = desiredTimes*fs;
indices = startAndEndIndices(1):startAndEndIndices(2);

% get the data from those indices and channels (ignoring ignored channels)
data = getiEEGData(dataName,channels,indices,pwfile);


%% Run spike detector
gdf = detectEKGSpikes(data.values,data.fs);
gdf = [zeros(length(gdf),1),gdf'];


if isempty(gdf) == 1
    fprintf('No spikes detected\n');
else
    fprintf('Detected %d spikes\n',size(gdf,1));
     % put it in seconds
    gdf(:,2) = gdf(:,2)/data.fs;
    
    % re-align time
    gdf(:,2) = gdf(:,2) + desiredTimes(1);
end


end