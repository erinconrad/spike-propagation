function distNearestSOZ(pt,cluster,whichPts)

% Save file location
[~,~,~,resultsFolder,~] = fileLocations;
destFolder = [resultsFolder,'clustering/stats/'];
mkdir(destFolder);


if isempty(whichPts) == 1
    for i = 1:length(pt)
        if isempty(pt(i).seq_matrix) == 0
            if size(cluster(i).bad_cluster) < cluster(i).k
                whichPts = [whichPts,i];
            end
        end
    end
    
    if isequal(whichPts,[1,4,6,7,8,9,12,14,15,16,17,18,19,20,22,24,25,27,30,31]) == 0
        error('Warning, not doing correct patients!\n');
    end
elseif whichPts == 100
    whichPts = [1,4,6,8,9,12,15,17,18,19,20,22,24,25,27,30,31];
end



for whichPt = whichPts
    fprintf('Doing %s\n',pt(whichPt).name);
    
    % Get patient parameters
    locs = pt(whichPt).electrodeData.locs(:,2:4);
    szTimes = pt(whichPt).newSzTimes;
    soz = pt(whichPt).newSOZChs;
    
    % Reorder seizure times if out of order
    oldSzTimes = szTimes;
    szTimes = sort(szTimes,1);
    if isequal(oldSzTimes,szTimes) == 0
        fprintf('WARNING!!! %s seizure times out of order\n',pt(whichPt).name);
    end
    
    % Combine nearly equal seizure times
    newIdx = 2;
    newSzTimes = [];
    newSzTimes(1,:) = szTimes(1,:);
    for j = 2:size(szTimes,1)
        if abs(szTimes(j,1)-szTimes(j-1,1)) < 10 && ...
                abs(szTimes(j,2)-szTimes(j-1,2))
           newIdx = newIdx - 1; 
           newSzTimes(newIdx,1) = min(szTimes(j,1),szTimes(j-1,1));
           newSzTimes(newIdx,2) = max(szTimes(j,2),szTimes(j-1,2));  
        else
           newSzTimes(newIdx,:) = szTimes(j,:);
        end
        newIdx = newIdx + 1;
    end
    
    if isequal(newSzTimes,szTimes) == 0
        fprintf('WARNING!!! %s had duplicate seizure times\n',pt(whichPt).name);
    end
    
    % Pull cluster info
    %{
    all_times_all = pt(whichPt).cluster.all_times_all; % all spike times
    all_spikes = pt(whichPt).cluster.all_spikes; % all spike channels
    all_locs = pt(whichPt).cluster.all_locs;
    k = pt(whichPt).cluster.k; % the number of clusters
    idx = pt(whichPt).cluster.idx; % the cluster index for every spike
    C = pt(whichPt).cluster.C; % the centroids of the clusters
    bad_cluster = pt(whichPt).cluster.bad_cluster; % which clusters are bad
    %}
    
    all_times_all = cluster(whichPt).all_times_all; % all spike times
    all_spikes = cluster(whichPt).all_spikes; % all spike channels
    all_locs = cluster(whichPt).all_locs;
    k = cluster(whichPt).k; % the number of clusters
    idx = cluster(whichPt).idx; % the cluster index for every spike
    C = cluster(whichPt).C; % the centroids of the clusters
    bad_cluster = cluster(whichPt).bad_cluster; % which clusters are bad
    
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
    C(bad_cluster,:) = [];
    
    
    % For each spike, get the distance between the spike and the nearest
    % SOZ
    soz = pt(whichPt).newSOZChs; 
    spike_dist = zeros(size(all_locs,1),1);
    for i = 1:size(all_locs,1)
        spike_dist(i) = min(vecnorm(all_locs(i,:) - locs(soz,:),2,2)); 
    end
    
    
    fprintf('Made it here\n');
    
    plot_times = all_times_all-min(all_times_all);
    scatter(plot_times/3600,spike_dist);
    hold on
    for j = 1:size(szTimes,1) 
        yl = ylim;
        plot([szTimes(j,1)-min(all_times_all) szTimes(j,1)-min(all_times_all)]/3600,...
            yl,'k','LineWidth',2);
    end
    
    
end



end