%% run_example_data

%{
The purpose of this script is to have a complete pipeline for running the
spike detector, spike sequence detector, and clustering algorithm on a set
of eeg data. This assumes you have a Matlab pt structure with eeg data. An
example structure can be downloaded along with this file.

This script then takes this eeg data and goes through each step in the
process, finally generating the spike clusters for viewing.

%}

%% Parameters
whichPt = 31; % This is the patient who has example data
surround_time = 3; % Seconds before and after the spike to plot


times = pt(whichPt).runTimes; % run times over which to run the spike detector
thresh = pt(whichPt).thresh; % threshold data for the spike detector
whichDetector = thresh.whichDetector; % which detector to use
example = 1; % a flag that we are using example data (otherwise this will try to find the data on ieeg, which requires an account)
fs = pt(whichPt).fs; % Sampling rate



%% Run spike detector
[gdf,~] = getSpikesSimple(pt,whichPt,times,whichDetector,thresh,0,example);
% gdf is an nx2 array, where there are n spikes, and the first column has
% the spike channel and the second column has the spike time. They are
% sorted by time.

pt(whichPt).gdf = gdf; % add the spike times to the pt structure

%% Plot an example spike
sp = 100; % change this to see other spikes

% Get spike time and ch
sp_time = gdf(sp,2);
sp_ch = gdf(sp,1);
elec_label = pt(whichPt).electrodeData.electrodes(sp_ch).name;
sp_idx = round((sp_time - pt(whichPt).eeg_data.times(1))*fs); % get the index of the spike
plot_indices = max(1,round(sp_idx - surround_time*fs)):...
    min(length(pt(whichPt).eeg_data.times),round(sp_idx + surround_time*fs)); % get the indices to plot
figure
set(gcf,'position',[110 452 987 352])
plot((plot_indices-plot_indices(1))/fs,pt(whichPt).eeg_data.values(plot_indices,sp_ch),'k') % plot the eeg data for the ch with the spike
hold on
scatter((sp_idx-plot_indices(1))/fs,pt(whichPt).eeg_data.values(sp_idx,sp_ch),80,'k','filled') % mark the spike time
yticklabels([])
xlim([0,surround_time*2])
xlabel('Time (s)')
title(sprintf('Spike time: %1.1f s, electrode %s',sp_time,elec_label))
set(gca,'fontsize',20)


%% Run spike sequence detector
% This detects multi-channel spike sequences
pt(whichPt).data = mainSequences(gdf,pt(whichPt).electrodeData, pt(whichPt).fs);

% This converts the sequences into an nch x nseq array, where nch is the
% number of channels and nseq is the number of sequences. This is mostly
% NaNs, but the non-NaN values are the spike times for the channels
% involved in each sequence
pt(whichPt).seq_matrix = ...
        makeSeqMatrix(pt(whichPt).data.sequences,length(pt(whichPt).channels),0);
    
%% Plot an example sequence
s = 5; % change this to see other sequences

seq = pt(whichPt).seq_matrix(:,s); % get which sequence
chs = find(~isnan(seq)); % get the spike channels in the sequence
times = seq(chs); % get spike times
[times,I] = sort(times); % sort by times
chs = chs(I);
sp_idx = (times(1) - pt(whichPt).eeg_data.times(1))*fs; % get time of first spike
plot_indices = max(1,round(sp_idx - surround_time*fs)):...
    min(length(pt(whichPt).eeg_data.times),round(sp_idx + surround_time*fs)); % center plot around the first spike

figure
set(gcf,'position',[110 452 987 352])
to_add = 0;
for i = 1:length(times)
    elec_label = pt(whichPt).electrodeData.electrodes(chs(i)).name;
    curr_sp_idx = round((times(i)- pt(whichPt).eeg_data.times(1))*fs); % get index of current spike
    plot((plot_indices-plot_indices(1))/fs,pt(whichPt).eeg_data.values(plot_indices,chs(i))+to_add,'k') % plot the data in that channel
    hold on
    scatter((curr_sp_idx-plot_indices(1))/fs,pt(whichPt).eeg_data.values(curr_sp_idx,chs(i))+to_add,'k') % show the spike
    text(6,to_add,sprintf('%s',elec_label),'fontsize',20)
    to_add = to_add - 2000; % to offset the lines so they don't blur into one another
end
yticklabels([])
xlim([0,surround_time*2])
xlabel('Time (s)')
title(sprintf('First spike time: %1.1f s',times(1)))
set(gca,'fontsize',20)
    
%% Get the optimal cluster number
% Note that this is imperfect because we are running the clustering
% algorithm over a short period of data. I tried to pick a section of data
% in which there seemed to be 3 separate spike clusters for example. When
% this is run, it will plot the sum-squared error (SSE) as a function of
% cluster number. Pick the cluster number the forms the inflection point of
% the graph to get the optimal cluster number. (It appears to be 3 for the
% example data).
getClusters(pt,whichPt,1,1,0,0); 

%% Run the clustering algorithm for K = 3 (looks like the inflection point on the SSE plot)
% This will run the clustering algorithm with a cluster number of 3. It
% will cluster the spikes based on their spike location. It will output a
% new structure "cluster" with the cluster identities of each spike. It
% will also plot 10 sample sequences for each cluster. These sequences can
% be manually reviewed to see if they appear to represent real spike
% sequences (or artifact that was incorrectly detected as a spike).

% One thing to note is that the definition of which is cluster 1, 2, 3,
% etc. will change each time. It is non-deterministic because of the nature
% of Matlab's k-means and it is NOT the case that the most populous cluster
% is cluster 1.
cluster = getClusters(pt,whichPt,0,1,1,0,3); 

%% Plot the cluster centroid locations
locs = pt(whichPt).electrodeData.locs(:,2:4); % Get all electrode locations
centroid_locs = cluster(whichPt).C; % Get centroid locations for each cluster
figure
colors = [0 0 1;1 0 0;0 1 0];
scatter3(locs(:,1),locs(:,2),locs(:,3),300,'k') % plot all electrodes
hold on
for i = 1:size(centroid_locs,1)
    scatter3(centroid_locs(i,1),centroid_locs(i,2),centroid_locs(i,3),400,...
        colors(i,:),'filled') % plot cluster centroids
end

%% Plot proportion of sequences in each cluster over time
% This defaults to not plot because it is likely not meaningful over such a
% short time period
if 1 == 1
    
    colors = [0 0 1;1 0 0;0 1 0]; %The colors (will work for up to 3 clusters)
    plot_times = cluster(whichPt).all_times_all; % plot times
    k = cluster(whichPt).k; % the number of clusters 
    idx = cluster(whichPt).idx; % the cluster index for every spike
    all_locs = cluster(whichPt).all_locs; % spike locations
    clusters = 1:k; % cluster numbers
    c_idx = zeros(size(idx,1),3); % initialize color index array

    % Get the times for spikes in each cluster
    for i = 1:length(clusters)
        clust{i} = plot_times(idx == clusters(i));
    end
    window = 30; % 30 second window over which to calculate moving average

    [sum_c,sum_times] = movingSumCounts(clust,plot_times,window); % get moving average
    totalSum = zeros(1,size(sum_times,2));
    
    % calculate moving average of total number of spikes
    for i = 1:length(clusters)
        totalSum = totalSum + sum_c(i,:);
    end
    
    % calculate proportion of spikes in each cluster at each time step
    prop_c = sum_c./totalSum;

    % Get colors for each cluster
    for i = 1:length(clusters)  
       c_idx(idx==clusters(i),:) = ...
           repmat(colors(i,:),sum(idx==clusters(i)),1);
    end

    % plot moving average of proportion of spikes in each cluster
    figure
    set(gcf,'position',[110 452 987 352])
    set(gcf,'position',[1 409 1440 320]);
    for i = 1:length(clusters)
        pl(i)= plot((sum_times-min(plot_times)),prop_c(i,:),...
            'color',colors((i),:),'LineWidth',2);
    hold on
    xlabel('Time (s)')
    ylabel(sprintf('Proportion in each cluster\n(moving average)'))
    legend({'Cluster 1','Cluster 2','Cluster 3'})
    set(gca,'fontsize',20)
    end
end
