function CLocations(pt,cluster,whichPts)

% Save file location
[~,~,~,resultsFolder,~] = fileLocations;
destFolder = [resultsFolder,'clustering/stats/'];
mkdir(destFolder);


if isempty(whichPts) == 1
    for i = 1:length(pt)
        if isempty(pt(i).seq_matrix) == 0
            whichPts = [whichPts,i];
        end
    end
end

for whichPt = whichPts
    
    fprintf('Doing %s\n',pt(whichPt).name);
    
     % Get patient parameters
    locs = pt(whichPt).electrodeData.locs(:,2:4);
    szTimes = pt(whichPt).newSzTimes;
    soz = pt(whichPt).newSOZChs;
    
    all_times_all = cluster(whichPt).all_times_all; % all spike times
    all_spikes = cluster(whichPt).all_spikes; % all spike channels
    all_locs = cluster(whichPt).all_locs;
    k = cluster(whichPt).k; % the number of clusters
    idx = cluster(whichPt).idx; % the cluster index for every spike
    C = cluster(whichPt).C; % the centroids of the clusters
    bad_cluster = cluster(whichPt).bad_cluster; % which clusters are bad
    allSpikes = cluster(whichPt).allSpikes; % am I taking all spikes???
    
    % Confirm that I do not have any ictal spikes
    t = find(any(all_times_all >= szTimes(:,1)' & all_times_all <= szTimes(:,2)',2));
    if isempty(t) == 0
        fprintf('WARNING: Remaining ictal spikes for %s!\n',pt(whichPt).name);
        all_times_all(t) = [];
        all_spikes(t) = [];
        all_locs(t,:) = [];
        idx(t) = [];
    end
    

    % Remove bad clusters
    bad_idx = find(ismember(idx,bad_cluster));
    all_times_all(bad_idx) = [];
    all_spikes(bad_idx) = [];
    all_locs(bad_idx,:) = [];
    idx(bad_idx) = [];
    clusters = 1:k; clusters(bad_cluster) = [];
    
    % Plot locations of spikes in main clusters
    for i = 1:length(clusters)
        
        % Get 100 random spikes in the cluster
        n_in_clust = sum(idx==i);
        all_in_clust = all_spikes(idx==i);
        spikes =  all_in_clust((randsample(n_in_clust,min(100,length(all_in_clust)))));
        spike_locs = locs(spikes,:);
        
        figure
        scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
        hold on
        scatter3(spike_locs(:,1),spike_locs(:,2),spike_locs(:,3),100,'r')
        pause
        close(gcf)
        
        
    end
    
    
end

end