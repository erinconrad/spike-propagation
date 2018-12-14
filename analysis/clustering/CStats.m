function stats = CStats(pt,cluster,whichPts)

%{ 

CStats
This is my cleaned up file for getting statistics on the cluster data 

%}

% Parameters
plotQI = 0;

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

allOutcome = [];
for whichPt = whichPts
    outcome(whichPt) = getOutcome(pt(whichPt).name);
    allOutcome = [allOutcome;outcome(whichPt)];
end

allCounts = [];
allPat = [];
allChunk = [];
chi_tables_plot = cell(max(whichPts),1);
p_plot = zeros(max(whichPts),1);


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
    
    %%
    %{
    ------
    Analysis 1: Does cluster distribution change from hour to hour? 
    ------
    %}
    test_t = 3600; % 60 minute chunks
    
    % Divide run into 60 minute chunks
    % This may result in some chunks that are empty because there was low
    % voltage data and so I skipped spike detection in this period, but
    % this should not affect the cluster distribution. 
    n_chunks = ceil((max(all_times_all) - min(all_times_all))/test_t);
    which_chunk = zeros(length(all_times_all),1);
    
    for i = 1:n_chunks
        
        % Get the time range for the chunk
        curr_times = [min(all_times_all) + (i-1)*test_t,...
           min(min(all_times_all) + i*test_t,max(all_times_all))];
       
        % Get the spike indices in that time chunk
        chunk_spikes = find(all_times_all >= curr_times(1) & ...
            all_times_all <= curr_times(2));
        
        % Label these spikes with a chunk index
        which_chunk(chunk_spikes) = i;
       
    end
    
    % Do a chi-squared test to see if the cluster distribution changes
    % across the 60 minute chunks
    [tbl_1,chi2_1,p_1,labels_1] = crosstab(which_chunk,idx);
    
    % Save information into patient struct
    stats(whichPt).cluster.hour.tbl = tbl_1;
    stats(whichPt).cluster.hour.chi2 = chi2_1;
    stats(whichPt).cluster.hour.p = p_1;
    stats(whichPt).cluster.hour.labels = labels_1;
    
    fprintf(['For %s, regarding whether 60 minute chunks\n have different cluster'...
    ' distributions,\n the p-value is %1.1e\n\n\n'],pt(whichPt).name,p_1);

    
    %%
    %{
    -------
    Analysis 2: Does cluster distribution differ between the pre-ictal and
    inter-ictal period?
    -------
    %}
    
    % Define important ranges
    preIcRange = [-60*60,-1*60]; % Between 1 hour to 1 minute before a seizure
    postIcTime = 60*60; % 60 minutes after a seizure
    % Interictal range will be anything else
    
    % Get times for QI purposes
    preIcTimesQI = [];
    interIcTimesQI = [];
    
    % Get all the pre-ictal spikes
    preIcClustIdx = [];
    
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
        
        preIcTimesQI = [preIcTimesQI;preIcTimes(1) preIcTimes(2)];
        
        % Get the indices of the spikes in this time range
        spike_idx = (all_times_all >= preIcTimes(1) & ...
            all_times_all <= preIcTimes(2));
        
        % Get the cluster indices
        preIcClustIdx = [preIcClustIdx;idx(spike_idx)];
        
    end
    
    % Get all the inter-ictal spikes
    interIcClustIdx = [];
    
    % Get times before the first seizure
    interIcTimes(1) = min(all_times_all) + postIcTime; % assume seizure right before the record started
    interIcTimes(2) = szTimes(1,1) + preIcRange(1) - 1; % Up to 60 minutes before first sz
    
    
    if interIcTimes(1) < interIcTimes(2)
        spike_idx = all_times_all >= interIcTimes(1) & ...
            all_times_all <= interIcTimes(2);
        interIcClustIdx = [interIcClustIdx; idx(spike_idx)];
        interIcTimesQI = [interIcTimesQI; interIcTimes(1) interIcTimes(2)];
    end
    
    % Loop through the seizures
    for j = 1:size(szTimes,1)-1
        interIcTimes(1) = szTimes(j,2) + postIcTime;
        interIcTimes(2) = szTimes(j+1,1) + preIcRange(1) - 1;
        
        if interIcTimes(1) >= interIcTimes(2), continue; end
        spike_idx = all_times_all >= interIcTimes(1) & ...
            all_times_all <= interIcTimes(2);
        interIcClustIdx = [interIcClustIdx; idx(spike_idx)]; 
        interIcTimesQI = [interIcTimesQI; interIcTimes(1) interIcTimes(2)];
    end
    
    % Get time after the last seizure
    interIcTimes(1) = szTimes(end,2) + postIcTime;
    interIcTimes(2) = max(all_times_all) + preIcRange(1) -1; % Assume seizure right after record ends
    
    if interIcTimes(1) < interIcTimes(2)
        spike_idx = all_times_all >= interIcTimes(1) & ...
            all_times_all <= interIcTimes(2);
        interIcClustIdx = [interIcClustIdx; idx(spike_idx)];
        interIcTimesQI = [interIcTimesQI; interIcTimes(1) interIcTimes(2)];
    end
    
    % Plot interical and preictal times next to seizure times for QI
    % purposes
    if plotQI == 1
        figure
        for j = 1:size(preIcTimesQI,1)
           area([preIcTimesQI(j,1) preIcTimesQI(j,2)]/3600,[1 1],'FaceColor','g');
           hold on
        end
        for j = 1:size(interIcTimesQI,1)
           area([interIcTimesQI(j,1) interIcTimesQI(j,2)]/3600,[1 1],'FaceColor','r');
           hold on
        end
        yl = ylim;
        for j = 1:size(szTimes,1)
           plot([szTimes(j,1) szTimes(j,1)]/3600,yl,'k--','LineWidth',5);
        end
    end
    
    % Do chi_2 to test if pre-ictal cluster distribution is different from
    % inter-ictal cluster distribution
    [tbl_2,chi2_2,p_2,labels_2] = crosstab([ones(size(preIcClustIdx));...
        2*ones(size(interIcClustIdx))],[preIcClustIdx;interIcClustIdx]);
    
    % Save information into patient struct
    stats(whichPt).cluster.preic.tbl = tbl_2;
    stats(whichPt).cluster.preic.chi2 = chi2_2;
    stats(whichPt).cluster.preic.p = p_2;
    stats(whichPt).cluster.preic.labels = labels_2;
    
    chi_tables_plot{whichPt} = tbl_2;
    p_plot(whichPt) = p_2;
    
    %{
    fprintf(['For %s, regarding whether the pre-ictal period\n has a different cluster'...
    ' distribution from the interictal period,\n the p-value is %1.1e\n\n'],pt(whichPt).name,p_2);
    %}
end

%% Fisher's method to combine p values
all_p = [];
for whichPt = whichPts
   all_p = [all_p;p_plot(whichPt)]; 
end

X_2 = -2 * sum(log(all_p));
sum_p = 1-chi2cdf(X_2,2*length(all_p));

fprintf('The group p value is %1.1e\n',all_p);

% double check
group_pval = fisher_pvalue_meta_analysis(all_p);

%% Bar graphs
figure
numcols = 10;
numrows = ceil(length(whichPts)/numcols);
[ha, pos] = tight_subplot(numrows,numcols,[.08 .01],[.05 .05],[.05 .01]); 
for j = 1:length(whichPts)
    axes(ha(j));
    
    % Get chi2 table
    tbl = chi_tables_plot{whichPts(j)};
    
    % Get proportions rather than absolute numbers
    prop = tbl./sum(tbl,2);
    
    b=bar(prop);
    xticklabels({'Pre-ic','Inter-ic'})
    
    % Plot p-value
    p = p_plot(whichPts(j));
    if p < 0.001
        textp = 'p < 0.001';
    else
        textp = sprintf('p = %1.3f',p);
    end
    hold on
    plot([1 2],[0.9 0.9]+0.01,'k')
    text(1.5,0.9+0.05,textp,'HorizontalAlignment','center',...
        'fontsize',15);
    
    legend_names = cell(size(tbl,2),1);
    for k = 1:length(legend_names)
       legend_names{k} = sprintf('Cluster %d',k);  
    end
    yticklabels([])
   % lgnd=legend(b,legend_names,'location','northwest');
   % set(lgnd,'color','none');
    title(sprintf('%s',pt(whichPts(j)).name),'fontsize',15);
    if j == 1
        ylabel('Proportion of sequences');
    end
    set(gca,'fontsize',15);
    
    
    
    
end

pause
print(gcf,[destFolder,'clustBar'],'-depsc');
eps2pdf([destFolder,'clustBar','.eps'])


%% Now correlate with outcome


save([destFolder,'stats.mat','stats']);


end