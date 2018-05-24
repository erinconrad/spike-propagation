function gdf = vMakeSegments(gdf,values,fs,vtime)

% This will take a gdf file which contains the time and location of
% detected spikes and it will output an array called delay, which is an nxm
% array where n is the number of spikes and m is the number of channels,
% and for each spike, the columns show the time in that segment at which
% the channel has its maximum; as well as an array called rms, which is the
% same size and for each spike, the collumns show the rms power for that
% channel in the spike

% calculate how many indices to go back and how many to go forward after
% the spike (vtime is a 1x2 array where the first number is a negative
% number representing how many seconds to go back and the second number is
% a positive number representing how many seconds to go forward)
nIdx = round(fs*vtime);

n_channels = size(values,2);


% get spike times
spikeTimes = gdf(:,2);

%% How to throw out overlapping spikes

% establish idxToRemove, a column of 1s for now
idxToKeep = ones(size(spikeTimes));

% the index of the last kept spike
origIdx = 1;

% loop through spikes
for i = 2:length(spikeTimes)
    
    % if the new spike is too close to the last kept spike
    if spikeTimes(i) - spikeTimes(origIdx) < (vtime(2)-vtime(1))
        
        % set idxToKeep of this index to 0, indicating we will toss it
        idxToKeep(i) = 0;
        
    else
    % if it's acceptable, then set this new spike as the last kept spike
        origIdx = i;
        
    end
end

% remove spike times that are too close to the last spike
spikeTimes(idxToKeep==0) = [];

% get the indices of the new spikes
spikeIdx = round(spikeTimes*fs);

%% Calculate delay and rms, the output variables

% Initialize the output variable delay, which for each spike segment gives
% the time relative to the initial index of each channel's maximum. 
delay = zeros(length(spikeIdx),n_channels);

% Initialize the output variable rms
rms = zeros(length(spikeIdx),n_channels);

% Loop through the spike indices
for i = 1:length(spikeIdx)
    
    % Loop through each channel
    for j = 1:n_channels
     
        % Get the values in that segment for that channel
        temp_signal = values(spikeIdx(i)+nIdx(1):spikeIdx(i)+nIdx(2),j);
        
        % Get the index corresponding to the maximum value the channel has
        % in that segment. This is the delay for that channel within that
        % segment
        [~,I] = max(temp_signal);
        delay(i,j) = I;
        
        % calculate root mean square for the channel
        rms(i,j) = sqrt(mean((temp_signal-mean(temp_signal)).^2));
    end
    
    % get minimum delay amongs the channels and subtract this from each
    % delay to make min 0
    delay(i,:) = delay(i,:) - min(delay(i,:));
    
    
end

gdf = struct;
gdf.delay = delay;
gdf.rms = rms;
gdf.spikeTimes = spikeTimes;

end