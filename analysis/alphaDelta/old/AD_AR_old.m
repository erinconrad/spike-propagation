function AD_AR_old(pt,cluster,power,whichPts)

%{
This calculates the correlation between the alpha delta power ratio and
features of spike sequences
%}

% Plot info about the model?
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

%% initialize
allP = [];
allT = [];
outcome_all = [];
temp_lobe_all = [];

all_b_soz = [];
all_p_soz = [];
all_b_SL = [];
all_p_SL = [];
all_t_soz = [];
all_t_SL = [];
all_b = [];
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
    
    %% Get most popular cluster
    popular = mode(idx);
    
   
    
    %% Get the sequence lengths of all spikes
    [seq_lengths,seq_times,seq_areas] = getSeqDist(pt,cluster,whichPt);
    
    bin_times = pt(whichPt).runTimes;
    prop_pop = zeros(size(bin_times,1),1);
    locs_bin = zeros(size(bin_times,1),3);
    soz_dist_bin = zeros(size(bin_times,1),1);
    num_spikes = zeros(size(bin_times,1),1);
    SD_bin = zeros(size(bin_times,1),1);
    std_z_bin = zeros(size(bin_times,1),1);
    prop_pop_chunk = zeros(size(bin_times,1),1);
    distNeeded = zeros(size(bin_times,1),1);
    SL_bin = zeros(size(bin_times,1),1); 
    SA_bin = zeros(size(bin_times,1),1); 
    n_s = zeros(size(bin_times,1),1);
    n_f = zeros(size(bin_times,1),1);
    
    % Run through bin times and get proportion of spikes in most popular
    % cluster for that bin.
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
        locs_bin(i,:) = mean(all_locs(whichSpikes,:),1);
        num_spikes(i) = length(whichSpikes);
        
        
        
        % Get mean sequence length in these times
        SL_bin(i) = mean(seq_lengths(seq_times > bin_times(i,1) & ...
            seq_times < bin_times(i,2)));
        
        % Get mean sequence area in these times
        SA_bin(i)  = mean(seq_areas(seq_times > bin_times(i,1) & ...
            seq_times < bin_times(i,2)));
        
        
        % Get mean distance from spike to nearest SOZ
        soz_dist = zeros(length(whichSpikes),1);
        
        % Get soz mean
        median_soz = median(soz,1);

        % Loop through all spikes in bin
        for j = 1:length(whichSpikes)

            % For each spike, get the distance between the spike and its
            % nearest SOZ
            soz_dist(j) = min(vecnorm(locs(all_spikes(whichSpikes(j)),:) - ...
                median_soz,2,2)); 
        end
        % Average that distance over all spikes in the bin
        soz_dist_bin(i) = mean(soz_dist);
        
        
        % Now get a measure of spatial dispersion for the bin
        n_spikes = length(whichSpikes); % number of spikes
        locs_sp = all_locs(whichSpikes,:); % location of every spike in the bin
        
        % Standard distance
        if size(locs_sp,1) < 20
            SD_bin(i) = nan;
        else
            SD_bin(i) = standardDistance(locs_sp);
        end
        std_z_bin(i) = std(locs_sp(:,3));
        
        % How many spikes are in most popular cluster for that chunk?
        prop_pop_chunk(i) = sum(whichClust == mode(whichClust))/length(whichClust);
        
        % distance from median to 60 percent of spikes in the chunk
        if size(locs_sp,1) < 10
            distNeeded(i) = nan;
        else
            distNeeded(i) = distToCaptureMost(locs_sp);
        end
        
    end
    
    all_ad = power(whichPt).alpha./power(whichPt).delta;
    mean_ad = nanmean(all_ad,1)';
    
    times = mean(bin_times,2);
    
    old_times = times;
    old_mean_ad = mean_ad;
    
    
    %% DO SOZ ANALYSIS
    fprintf('Doing SOZ for %s\n',pt(whichPt).name);
    %% Remove nans for distance from SOZ analysis
    nan_times = find(isnan(soz_dist_bin));
    times = old_times;
    mean_ad = old_mean_ad;
    
    times(nan_times) = [];
    mean_ad(nan_times) = [];
    soz_dist_bin(nan_times) = [];
    
    % Do model
    Y = soz_dist_bin;
    
    % X is the predictor. The first component of X is the alpha-delta
    % ratio, the predictor I am interested in. The next component is a
    % constant error term. The third component is just what time it is,
    % reflecting a linear trend with time. The last is the categorical
    % variable representing the hour of the day, reflecting a cyclical Q24
    % hour trend.
    %X = [mean_ad ones(size(mean_ad)) times cat_hours];
    X = [mean_ad ones(size(mean_ad))];
    
    %[p,t,b] = AR_model(X,Y,plotInfo);
    [p,t,b] = determine_order(X,Y,plotInfo);
    
    all_b_soz = [all_b_soz;b];
    all_p_soz = [all_p_soz;p];
    all_t_soz = [all_t_soz;t];
    
    
    
    %% Now do SL
    fprintf('Doing SL for %s\n',pt(whichPt).name);
    
    nan_times = find(isnan(SL_bin));
    times = old_times;
    mean_ad = old_mean_ad;
    
    times(nan_times) = [];
    mean_ad(nan_times) = [];
    SL_bin(nan_times) = [];
    
    %% Output x and y into a csv to be loaded into R to do the glarma model
    M = [mean_ad,SL_bin,ones(size(SL_bin))];
    fname = [resultsFolder,'for_r/',pt(whichPt).name,'_SL.csv'];
    csvwrite(fname,M)
    
    % Do model
    
    Y = SL_bin;
    
    % X is the predictor. The first component of X is the alpha-delta
    % ratio, the predictor I am interested in. The next component is a
    % constant error term. The third component is just what time it is,
    % reflecting a linear trend with time. The last is the categorical
    % variable representing the hour of the day, reflecting a cyclical Q24
    % hour trend.
    %X = [mean_ad ones(size(mean_ad)) times cat_hours];
    X = [mean_ad ones(size(mean_ad))];
    
    %[p,t,b] = AR_model(X,Y,plotInfo);
    [p,t,b] = determine_order(X,Y,plotInfo);
    
    all_b_SL = [all_b_SL;b];
    all_p_SL = [all_p_SL;p];
    all_t_SL = [all_t_SL;t];
   
    
    
    %% DO PROP-POP ANALYSIS
    fprintf('Doing prop-pop for %s\n',pt(whichPt).name);
    
    % Skip if only one cluster
    
    if length(clusters) == 1
        fprintf('One cluster for %s, skipping\n',pt(whichPt).name);
        p = 1;
        allP =[allP;p];
        allT = [allT;nan];
        all_b = [all_b;nan];
        all_b_glm = [all_b_glm;nan];
        all_p_glm = [all_p_glm;nan];
        continue
    end
    
    % Remove nan time points for prop-pop analysis
    nan_times = find(isnan(prop_pop));
    times = old_times;
    mean_ad = old_mean_ad;
    
    
    times(nan_times) = [];
    mean_ad(nan_times) = [];
    prop_pop(nan_times) = [];
    n_s(nan_times) = [];
    n_f(nan_times) = [];
    
    
    %{
    % Get hours
    hours = floor(mod(times,24*3600)/3600) + 1;
    
    % Get hours as a categorical variable
    cat_hours = dummyvar(hours);
    %}
    
    % Y is the response variable, the proportion of sequences in the most
    % popular cluster
    Y = prop_pop;
    
    % X is the predictor. The first component of X is the alpha-delta
    % ratio, the predictor I am interested in. The next component is a
    % constant error term. The third component is just what time it is,
    % reflecting a linear trend with time. The last is the categorical
    % variable representing the hour of the day, reflecting a cyclical Q24
    % hour trend.
    %X = [mean_ad ones(size(mean_ad)) times cat_hours];
    X = [mean_ad ones(size(mean_ad))];
   
    %% Output x and y into a csv to be loaded into R to do the glarma model
    M = [mean_ad,n_s,n_f,ones(size(prop_pop)),prop_pop];
    fname = [resultsFolder,'for_r/',pt(whichPt).name,'.csv'];
    csvwrite(fname,M)
    
    % Do model
    [p,t,b] = determine_order(X,Y,plotInfo);
    %[p,t,~] = AR_model(X,Y,plotInfo);
    
    allP = [allP;p];
    allT = [allT;t];
    all_b = [all_b;b];
    
    %% GLM model
    [b_glm,~,stats_glm] = glmfit(mean_ad,[n_s n_s+n_f],'binomial', 'link', 'logit');
    yfit = glmval(b_glm, mean_ad, 'logit', 'size', n_s+n_f);
    p_glm = stats_glm.p(2);
    all_p_glm = [all_p_glm;p_glm];
    all_b_glm = [all_b_glm;b_glm(2)];
    
    
    %% Plot partial correlation of residuals
    %{
    parcorr(stats_glm.resid)
    pause
    close(gcf)
    %}
    
    %% Linear model
    [b_lin,~,~,~,stats_lin] = regress(n_s,X);
    yfit_lin = X*b_lin;
    
    
    
    %% Plot comparison
   %{
    figure
    %plot(mean_ad,n_s./(n_s+n_f),'o',mean_ad,yfit./(n_s+n_f),'-');
    plot(n_s./(n_s+n_f),'o')
    hold on
    plot(yfit./(n_s+n_f))
    pause
    close(gcf);
   %}
    
   %{
    figure
    set(gcf,'position',[100 100 1000 500])
    subplot(2,1,1)
    plot(n_s./(n_s+n_f),'o')
    hold on
    plot(yfit./(n_s+n_f))
    legend('Real proportion in most popular cluster','Logistic')
    
    subplot(2,1,2)
    plot(n_s./(n_s+n_f),'o')
    hold on
    plot(yfit_lin./(n_s+n_f))
    legend('Real proportion in most popular cluster','Linear')
    pause
    close(gcf)
    %}
    
end


X_2 = -2 * sum(log(allP));
sum_p = 1-chi2cdf(X_2,2*length(allP));
fprintf(['There are %d with significant correlation.\nThe combined p-value'...
    'is %1.1e.\n'],sum(allP<0.05/length(allP)),sum_p);

%% Make table of p-values and t-statistics
p_text = getPText(allP);
allT_text = num2str(allT,3);
all_b_text = num2str(all_b,3);
table(char(names),char(all_b_text),allT_text,char(p_text))

all_b_glm_text = num2str(all_b_glm,3);
all_p_glm_text = num2str(all_p_glm,3);
table(char(names),char(all_b_glm_text),char(all_p_glm_text))

%% Test that t significantly different from zero for prop-pop
[~,p,ci,stats] = ttest(allT);
fprintf(['P value for consistent change in proportion of spikes in\n'...
    'most popular cluster is %1.1e, tstat, %1.2f.\n'],p,stats.tstat);

%% Test that t significantly different from zero for SOZ
[~,p,ci,stats] = ttest(all_t_soz);
fprintf('P value for SOZ is %1.1e, tstat, %1.2f.\n',p,stats.tstat);
changePos = find(allP<0.05/length(allP));
all_t_soz_changePos = all_t_soz(changePos);
all_p_soz_changePos = all_p_soz(changePos);
table(names(changePos),all_t_soz_changePos,all_p_soz_changePos)

%% Correlation between change in location and outcome
%
changeLoc = allP < 0.05/length(allP);
[p,info] = correlateClinically(changeLoc,outcome_all,'bin','num',0);
fprintf(['The p-value for Wilcoxon rank sum comparing outcome between\n'...
    'patients with change and those without is:\n'...
    'p = %1.2e\nranksum = %1.1f\n'],p,info.stats.ranksum);


%% Correlation between change in location and temp lobe
[p,info] = correlateClinically(changeLoc,temp_lobe_all,'bin','bin',0);
fprintf(['The p-value for chi squared comparing temporal vs non-temporal lobe between\n'...
    'patients with change and those without is:\n'...
    'p = %1.2e\nchi2 =  %1.2f\n'],p,info.chi2);


%% Table of stuff
p_sleep_t = getPText(allP);
T = table(p_sleep_t);


%% Test that t-stat significantly different from zero for SL
[~,p,ci,stats] = ttest(all_t_SL);
fprintf('P value for SL is %1.1e, tstat %1.2f.\n',p,stats.tstat);
changePos = find(allP<0.05/length(allP));
all_t_SL_changePos = all_t_SL(changePos);
all_p_SL_changePos = all_p_SL(changePos);
table(names(changePos),all_t_SL_changePos,all_p_SL_changePos)

%% P

end


function p_text_cell = getPText(p_array)

    p_text_cell = cell(length(p_array),1);
    for i = 1:length(p_array)
        if p_array(i) < 0.001
            p_text_cell{i} = '<0.001';
        else
            p_text_cell{i} = sprintf('%1.3f',p_array(i));
        end
    
        if p_array(i) < 0.001/length(p_array)
            p_text_cell{i} = [p_text_cell{i},'***'];
        elseif p_array(i) < 0.01/length(p_array)
            p_text_cell{i} = [p_text_cell{i},'**'];
        elseif p_array(i) < 0.05/length(p_array)
            p_text_cell{i} = [p_text_cell{i},'*'];
        end
    
    end
        
        
end

