function glm_model(pt,cluster,power,whichPts)


%{
This script determines whether there is a relationship between the
proportion of spikes in the predominant cluster and the alpha delta ratio,
assuming a GLM model for the relationship, rather than a GLARMA model (for
which you would need to use R).

%}

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


%% initialize
names = {};
all_b_glm = [];
all_p_glm = [];

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
    
    
    %% Initialize things
    bin_times = pt(whichPt).runTimes;
    n_s = zeros(size(bin_times,1),1);
    n_f = zeros(size(bin_times,1),1);
    prop_pop = zeros(size(bin_times,1),1);
    
    % Get most popular cluster
    popular = mode(idx);
    
    %% Run through bin times and get proportion of spikes in most popular cluster for that bin.
    for i = 1:size(bin_times,1)
        
        
        
        
        % Cluster identities of spikes in between those times
        whichClust = idx(all_times_all > bin_times(i,1) & ...
            all_times_all < bin_times(i,2));
        
        prop_pop(i) = sum(whichClust == popular)/length(whichClust);
        n_s(i) = sum(whichClust == popular);
        n_f(i) = sum(whichClust ~= popular);

        
        % Get spike locs
        whichSpikes = find(all_times_all > bin_times(i,1) & ...
            all_times_all < bin_times(i,2));
        
    end
    
    %% Get alpha delta ratio
    all_ad = power(whichPt).alpha./power(whichPt).delta;
    mean_ad = nanmean(all_ad,1)';
    
    times = mean(bin_times,2);
    
     %% DO PROP-POP ANALYSIS
    fprintf('Doing prop-pop for %s\n',pt(whichPt).name);
    
    
    % Skip if only one cluster
    if length(clusters) == 1
        fprintf('One cluster for %s, skipping\n',pt(whichPt).name);
        all_b_glm = [all_b_glm;nan];
        all_p_glm = [all_p_glm;nan];
        continue
    end
    
    % Remove nan time points for prop-pop analysis
    nan_times = find(isnan(prop_pop));
    
    times(nan_times) = [];
    mean_ad(nan_times) = [];
    prop_pop(nan_times) = [];
    n_s(nan_times) = [];
    n_f(nan_times) = [];
    
    % GLM model
    [b_glm,~,stats_glm] = glmfit(mean_ad,[n_s (n_s+n_f)],'binomial', 'link', 'logit');
    yfit = glmval(b_glm, mean_ad, 'logit', 'size', (n_s+n_f));
    p_glm = stats_glm.p(2);
    all_p_glm = [all_p_glm;p_glm];
    all_b_glm = [all_b_glm;b_glm(2)];
    
end

all_b_glm_text = num2str(all_b_glm,3);
all_p_glm_text = num2str(all_p_glm,3);

if length(whichPts) > 1
    table(char(names),char(all_b_glm_text),char(all_p_glm_text),'VariableNames',...
    {'Patient','Model_coefficient','p_value'})
else
    fprintf('Model coefficient correlating proportion of spikes in\npredominant cluster with alpha delta ratio is\n%s and p value is %1.3e\n',char(all_b_glm_text),all_p_glm);
end

end