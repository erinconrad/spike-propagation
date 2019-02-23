function rate_AD(pt,cluster,power,whichPts)



plotInfo = 0;

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
p1 = genpath(scriptFolder);
addpath(p1);
destFolder = [resultsFolder,'alphaDelta/plots/'];
mkdir(destFolder);

%% Get which patients
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

all_b = [];
all_p = [];
all_t = [];
names = {};
outcome_all = [];
temp_lobe_all = [];

%% Loop through patients
for whichPt = whichPts
    
    fprintf('Doing %s\n',pt(whichPt).name);
    %look = mean(power(whichPt).ad_rat,2);
    
    names = [names;pt(whichPt).name];
    
    saveFolder = [destFolder,pt(whichPt).name,'/'];
    mkdir(saveFolder)
    fprintf('Doing %s\n',pt(whichPt).name);
    szTimes = pt(whichPt).newSzTimes;
    locs = pt(whichPt).electrodeData.locs(:,2:4);
    soz = locs(pt(whichPt).newSOZChs,:);
    
    %% outcomes
    outcome = getOutcome(pt,whichPt);
    outcome_all = [outcome_all,outcome];
    
    %% SOZ
    szOnsetText = pt(whichPt).clinical.seizureOnset;
    if contains(szOnsetText,'TL') == 1
        tempLobe = 1;
    else
        tempLobe = 0;
    end
    temp_lobe_all = [temp_lobe_all,tempLobe];

    
    %% Reorder seizure times if out of order
    oldSzTimes = szTimes;
    szTimes = sort(szTimes,1);
    if isequal(oldSzTimes,szTimes) == 0
        fprintf('WARNING!!! %s seizure times out of order\n',pt(whichPt).name);
    end
    
    %% Combine nearly equal seizure times
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
    
    %% Get clustering info
    all_times_all = cluster(whichPt).all_times_all; % all spike times
    all_spikes = cluster(whichPt).all_spikes; % all spike channels
    all_locs = cluster(whichPt).all_locs;
    k = cluster(whichPt).k; % the number of clusters
    idx = cluster(whichPt).idx; % the cluster index for every spike
    C = cluster(whichPt).C; % the centroids of the clusters
    bad_cluster = cluster(whichPt).bad_cluster; % which clusters are bad
    
    %% Confirm that I do not have any ictal spikes
    t = find(any(all_times_all >= szTimes(:,1)' & all_times_all <= szTimes(:,2)',2));
    if isempty(t) == 0
        fprintf('WARNING: Remaining ictal spikes for %s!\n',pt(whichPt).name);
        all_times_all(t) = [];
        all_spikes(t) = [];
        all_locs(t,:) = [];
        idx(t) = [];
    end
    

    %% Remove bad clusters
    bad_idx = find(ismember(idx,bad_cluster));
    all_times_all(bad_idx) = [];
    all_spikes(bad_idx) = [];
    all_locs(bad_idx,:) = [];
    idx(bad_idx) = [];
    clusters = 1:k; clusters(bad_cluster) = [];
    
    bin_times = pt(whichPt).runTimes;
    counts = zeros(size(bin_times,1),1);
    
    % Run through bin times and get proportion of spikes in most popular
    % cluster for that bin.
    for i = 1:size(bin_times,1)
        counts(i) = sum(all_times_all > bin_times(i,1) & ...
            all_times_all < bin_times(i,2)); 
    end
    
    all_ad = power(whichPt).alpha./power(whichPt).delta;
    mean_ad = nanmean(all_ad,1)';
    
    times = mean(bin_times,2);
    old_times = times;
    old_mean_ad = mean_ad;
    
    nan_times = find(isnan(mean_ad));
    counts(nan_times) = [];
    mean_ad(nan_times) = [];
    
    
    %% Do model
    Y = counts;
    X = [mean_ad ones(size(mean_ad))];
    [p,t,b] = determine_order(X,Y,plotInfo);
    
    all_b = [all_b;b];
    all_p = [all_p;p];
    all_t = [all_t;t];
    
    %{
     % Do model
    Y = counts;
    
    % X is the predictor. The first component of X is the alpha-delta
    % ratio, the predictor I am interested in. The next component is a
    % constant error term. The third component is just what time it is,
    % reflecting a linear trend with time. The last is the categorical
    % variable representing the hour of the day, reflecting a cyclical Q24
    % hour trend.
    %X = [mean_ad ones(size(mean_ad)) times cat_hours];
    X = [mean_ad ones(size(mean_ad))];
    
    
    [p,t,b] = AR_model(X,Y,plotInfo);
    
    all_b = [all_b;b];
    all_p = [all_p;p];
    all_t = [all_t;t];
    
    %}
end

[h,p,~,stats] = ttest(all_t)

all_b
all_p



table(char(names(all_p<0.0025)),all_b(all_p<0.0025),all_p(all_p<0.0025))

end