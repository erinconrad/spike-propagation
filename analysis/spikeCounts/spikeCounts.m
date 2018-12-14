function spikeCounts(pt,whichPts)

% Parameters
doPlots = 1;

[~,~,scriptFolder,resultsFolder,~] = fileLocations;
p1 = genpath(scriptFolder);
addpath(p1);

destFolder = [resultsFolder,'spikeCounts/'];
mkdir(destFolder)

if isempty(whichPts) == 1
    for i = 1:length(pt)
        if isempty(pt(i).seq_matrix) == 0
            whichPts = [whichPts,i];
        end
    end
end

for whichPt = whichPts
    
    %% Patient parameters
    fprintf('Doing %s\n',pt(whichPt).name);
    szTimes = pt(whichPt).newSzTimes;
    saveFolder = [destFolder,pt(whichPt).name,'/'];
    mkdir(saveFolder);
    
    %% Get all sequences
    seq_matrix = pt(whichPt).seq_matrix;
    
    %% Remove ictal sequences
    first_time = min(seq_matrix,[],1);
    t = (any(first_time >= (szTimes(:,1)-repmat(60,size(szTimes,1),1)) ...
        & first_time <= szTimes(:,2),2));
    seq_matrix(:,t) = [];
    fprintf('Removed %d ictal spikes \n',sum(t));
    
    
    %% Get spikes
    spike_times = [];
    spike_chs = [];
    seq_index = [];
    for i = 1:size(seq_matrix,2)
        nonan = find(~isnan(seq_matrix(:,i)));
        times = seq_matrix(nonan,i);
        
        % Resort by time
        [~,I] = sort(times);
        nonan = nonan(I);
        
        spike_times = [spike_times;seq_matrix(nonan,i)];
        spike_chs = [spike_chs;nonan];
        seq_index = [seq_index;i*ones(length(nonan),1)];
    end
    
    %% Remove chunks of missing data
    
    time_chunks = getUnremovedTimes(pt,whichPt);

    
    %{
    ------
    Analysis 1: Does spike count change from hour to hour? 
    ------
    %}
    
    test_t = 3600; % 60 minute chunks
    
    % Get total time over all unremoved time chunks
    total_time = sum(diff(time_chunks,1,2));
    n_bins = ceil(total_time/test_t);
    num_spikes = zeros(nbins,1);
    
    % Divide the time_chunks into n_nbins
    start_time = time_chunk(1,1);
    start_chunk = 1;
    curr_bin_time = 0;
    % Loop through bins
    for i = 1:n_bins
        
        % Loop through chunks
        for j = start_chunk:size(time_chunks,1)
            
            % Start with the start time, and then go up to whatever is
            % smallest - the end of the chunk or the start time + remaining
            % bin time
            curr_times = [start_time ...
                min(start_time + (test_t-curr_bin_time),time_chunks(j,2))];
            
            % Add the spikes that occur during that time
            num_spikes(i) = num_spikes(i) +...
                sum(spike_times >= curr_times(1) & ...
                spike_times <= curr_times(2));
            
            % Increase how much bin time we've reached
            curr_bin_time = curr_bin_time + curr_times(2) - curr_times(1);
            
            % Redefine start time and chunk
            start_time = curr_times(2);
            start_chunk = j;
            
            if curr_bin_time > test_t
                break
            end
            
        end
        
    end
    
    % Now, as a result of this, I have n_bins that have 3600 seconds of NON
    % MISSING DATA, and so I can compare the spike counts across these bins
    % to see if there is a change over time
    
    % the expected number of counts per bin is the total number of counts
    % divided by the total number of bins
    expected = sum(num_spikes)/n_bins;
    chi_sq = sum((num_spikes - expected).^2/expected);
    df = n_bins-1;
    p = 1-chi2cdf(chi_sq,df);
    
    %{
    ------
    Analysis 1: Does spike count differ between pre-ictal and interictal? 
    ------
    %}
    
    % NEED TO THINK ABOUT WHETHER THIS MAKES SENSE
    
    % Here I think my approach will be to get the pre-ictal and inter-ictal
    % spikes and then to throw out any that are in the removed times
    
    % Define important ranges
    preIcRange = [-60*60,-1*60]; % Between 1 hour to 1 minute before a seizure
    postIcTime = 60*60; % 60 minutes after a seizure
    % Interictal range will be anything else
   
    
    % Get all the pre-ictal spikes
    preIcIdx = [];
    preTimes = [];
    
    % Loop through seizures
    for j = 1:size(szTimes,1)
        
        % Get range of pre-ictal times
        preIcTimes = szTimes(j,1) + preIcRange;
        
        % Shorten it if it is in the post-ictal period for the prior
        % seizure
        if j > 1
            preIcTimes(1) = max(preIcTimes(1),szTimes(j-1,2) + postIcTime);
        end
        
        % Shorten it if too close to the beginning of the record
        preIcTimes(1) = max(preIcTimes(1),min(all_times_all));
        
        % Skip the pre-ictal period if now there are no times left
        if preIcTimes(1) >= preIcTimes(2), continue; end
        
        preTimes = [preTimes; preIcTimes(1) preIcTimes(2)];
  
        % Get the indices of the spikes in this time range
        spike_idx = (spike_times >= preIcTimes(1) & ...
            spike_times <= preIcTimes(2));
        
        % Get the indices
        preIcIdx = [preIcIdx;(spike_idx)];
        
    end
    
    % Now remove those that aren't in the time_chunks array
    preIcIdx = preIcIdx(any(spike_times(preIcIdx)' > time_chunks(:,1) & ...
        spike_times(preIcIdx)' < time_chunks(:,2)));
    
    % Break up preTimes
    preTimesT = preTimes';
    time_chunksT = time_chunks';
    
    unremovedPreTimes = range_intersection(preTimesT,time_chunksT);
    
    % Get all the inter-ictal spikes
    interIcIdx = [];
    interTimes = [];
    
    % Get times before the first seizure
    interIcTimes(1) = min(spike_times) + postIcTime; % assume seizure right before the record started
    interIcTimes(2) = szTimes(1,1) + preIcRange(1) - 1; % Up to 60 minutes before first sz
    
    
    if interIcTimes(1) < interIcTimes(2)
        spike_idx = spike_times >= interIcTimes(1) & ...
            spike_times <= interIcTimes(2);
        interIcIdx = [interIcIdx; (spike_idx)];
        interTimes = [interTimes; interIcTimes(1) interIcTimes(2)];
    end
    
    % Loop through the seizures
    for j = 1:size(szTimes,1)-1
        interIcTimes(1) = szTimes(j,2) + postIcTime;
        interIcTimes(2) = szTimes(j+1,1) + preIcRange(1) - 1;
        
        if interIcTimes(1) >= interIcTimes(2), continue; end
        spike_idx = spike_times >= interIcTimes(1) & ...
            spike_times <= interIcTimes(2);
        interIcIdx = [interIcIdx; idx(spike_idx)]; 
        interTimes = [interTimes; interIcTimes(1) interIcTimes(2)];
    end
    
    % Get time after the last seizure
    interIcTimes(1) = szTimes(end,2) + postIcTime;
    interIcTimes(2) = max(spike_times) + preIcRange(1) -1; % Assume seizure right after record ends
    
    if interIcTimes(1) < interIcTimes(2)
        spike_idx = spike_times >= interIcTimes(1) & ...
            spike_times <= interIcTimes(2);
        interIcIdx = [interIcIdx; (spike_idx)];
        interTimes = [interTimes; interIcTimes(1) interIcTimes(2)];
    end
    

    % Now remove those that aren't in the time_chunks array
    interIcIdx = interIcIdx(any(spike_times(interIcIdx)' > time_chunks(:,1) & ...
        spike_times(interIcIdx)' < time_chunks(:,2)));
    
    % Now remove times
    
    
end


end