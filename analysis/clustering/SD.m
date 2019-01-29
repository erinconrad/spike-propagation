function SD(pt,cluster,whichPts)

%% Define which patients I am doing

if isempty(whichPts) == 1
    for i = 1:length(pt)
        if isempty(pt(i).seq_matrix) == 0
            if size(cluster(i).bad_cluster) < cluster(i).k
                whichPts = [whichPts,i];
            end
        end
    end
    
    % I should be doing all 20 of these patients
    if isequal(whichPts,[1,4,6,7,8,9,12,14,15,16,17,18,19,20,22,24,25,27,30,31]) == 0
        error('Warning, not doing correct patients!\n');
    end
end

count = 0;
names = cell(length(whichPts),1);
SD_all = zeros(length(whichPts),1);
outcome = zeros(length(whichPts),1);
TL = zeros(length(whichPts),1);
SD_elecs_all = zeros(length(whichPts),1);


for whichPt = whichPts
    
    count = count + 1;
    
    fprintf('Doing %s\n',pt(whichPt).name);
    
    names{count} = pt(whichPt).name;
    outcome(count) = getOutcome(pt,whichPt);
    onset = pt(whichPt).clinical.seizureOnset;
    
    if contains(onset,'TL') == 1
        TL(count) = 1;
    else
        TL(count) = 0;
    end
    
    %% Get patient parameters
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
    
    %% Get cluster info
    
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
    
    %% Get SD over all time
    SD = standardDistance(all_locs);
    SD_all(count) = SD;
    
    SD_elecs_all(count) = standardDistance(locs);
    
    %% Divide run into hours and get SD per hour
    window = 3600;
    n_bins = ceil((all_times_all(end)-all_times_all(1))/window);
    SD_hour = zeros(n_bins,1);
    times_hour = zeros(n_bins,1);
    mean_loc = zeros(n_bins,1);
    range_bin = zeros(n_bins,1);
    for i = 1:n_bins
        t_range = [(i-1)*window + all_times_all(1), ...
            min(all_times_all(end),i*window + all_times_all(1))];
        which_spikes = find(all_times_all >= t_range(1) & ...
            all_times_all < t_range(2));
        which_locs = all_locs(which_spikes,:);
        if isempty(which_locs) == 1
            SD_hour(i) = nan;
            mean_loc(i) = nan;
            range_bin(i) = nan;
        else
            %SD_hour(i) = standardDistance(which_locs);
            SD_hour(i) = sqrt(sum((which_locs(:,3)-...
                mean(which_locs(:,3))).^2))/...
                length(which_locs);
            mean_loc(i) = mean(which_locs(:,3));
            range_bin(i) = max(which_locs(:,3)) - min(which_locs(:,3));
        end
        
        times_hour(i) = (t_range(1) + t_range(2))/2;
    end
    
    %% Do plot
    figure
    subplot(3,1,1)
    plot(times_hour/3600,SD_hour);
    
    subplot(3,1,2)
    scatter(all_times_all/3600,all_locs(:,3))
    
    subplot(3,1,3)
    plot(times_hour/3600,range_bin)
    
end

[p1,info1] = correlateClinically(SD_all,outcome,'num','num',0);

[p2,info2] = correlateClinically(SD_all,TL,'num','bin',0);



end