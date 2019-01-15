function stats = CStats(pt,cluster,whichPts)

%{ 

CStats
This is my cleaned up file for getting statistics on the cluster data 

%}

% Parameters
plotQI = 1;
intericTime = 1;
doPermute = 0;

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
elseif whichPts == 100
    whichPts = [1,4,6,8,9,12,15,17,18,19,20,22,24,25,27,30,31];
end

allCounts = [];
allPat = [];
allChunk = [];
chi_tables_plot = cell(max(whichPts),1);
p_plot = zeros(max(whichPts),1);
p_change_time = zeros(max(whichPts),1);

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
    
    p_change_time(whichPt) = p_1;
    
    %%
    %{
    -------
    Analysis 2: Does cluster distribution differ between the pre-ictal and
    inter-ictal period?
    -------
    %}
    
    % Define important ranges
    preIcRange = [-60*60,-1*60]; % Between 1 hour to 1 minute before a seizure
    postIcTime = 60*60*intericTime; % 60 minutes after a seizure versus 4 hours
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
        pause
        close(gcf)
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
    
    if doPermute == 1
    
        %{
        Permutation test approach: if we randomly swap interictal and pre-ictal
        identities, how often do we get as significant a p-value? 

        My concern with this is the temporal dependence... if there is an hour
        to hour change and I am grouping things by hour, then this may be
        feeding the difference...

        %}
        n_preIc = size(preIcClustIdx,1);
        n_interIc = size(interIcClustIdx,1);
        nboot = 1e3;
        chi2_boot = zeros(nboot,1);
        for ib = 1:nboot
            if mod(ib,100) == 0
                fprintf('Doing %d of %d\n',ib,nboot);
            end
            preIcVsInterIc = ones(n_preIc+n_interIc,1);

            % Get a random sample, equal to the true number of interictal
            % spikes, and define those to be interictal
            fakeInterIc = randsample(n_preIc+n_interIc,n_interIc);
            preIcVsInterIc(fakeInterIc) = 2;

            [~,chi2_boot(ib),~,~] = crosstab(preIcVsInterIc,[preIcClustIdx;interIcClustIdx]);
        end
        
        s_chi = sort(chi2_boot);
        diff_s_chi = s_chi-chi2_2;
        allLarger = find(diff_s_chi>0);
        if isempty(allLarger) == 1
            p_boot = 0;
        else
            firstLarger = allLarger(1);
            p_boot = (nboot-firstLarger)/nboot;
        end
        
        fprintf('By bootstrap, pre-ictal vs inter-ictal cluster distribution diff is:\n%1.1e\n',...
            p_boot);

        figure
        scatter(1:nboot,sort(chi2_boot))
        hold on
        plot(xlim,[chi2_2 chi2_2])
        pause
        close(gcf)

        
    end
end

%% Fisher's method to combine p values
all_p = [];
for whichPt = whichPts
   all_p = [all_p;p_plot(whichPt)]; 
end

X_2 = -2 * sum(log(all_p));
sum_p = 1-chi2cdf(X_2,2*length(all_p));

fprintf('The group p value is %1.1e\n',sum_p);

% double check
%group_pval = fisher_pvalue_meta_analysis(all_p);

%% Get max number of clusters amongst patients (for making legend)
ptMax = 9;
maxNum = 1;
for whichPt = whichPts
    n_cl = size(chi_tables_plot{whichPt},2);
    if n_cl > maxNum
        maxNum = n_cl;
        ptMax = whichPt;
    end
end

%% Bar graphs
figure
set(gcf,'Position',[99 26 1264 686]);
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
        flag = '***';
    elseif p < 0.01
        flag = '**';
    elseif p < 0.05
        flag = '*';
    else
        flag='';
    end
    if p < 0.001
        textp = sprintf('p < 0.001%s',flag);
    else
        textp = sprintf('p = %1.3f%s',p,flag);
    end
    hold on
    plot([1 2],[0.9 0.9]+0.01,'k')
    text(1.5,0.9+0.05,textp,'HorizontalAlignment','center',...
        'fontsize',15);
    
    legend_names = cell(size(tbl,2),1);
    for k = 1:length(legend_names)
       legend_names{k} = sprintf('Cluster %d',k);  
    end
    if mod(j,10) ~= 1
        yticklabels([])
    end
   
    if whichPts(j) == ptMax
       lgnd=legend(b,legend_names,'Position',...
           [pos{length(whichPts)+1}(1)+0.023 pos{length(whichPts)+1}(2)+0.31,...
           0.03 0.1],...
           'units','normalized');
       set(lgnd,'color','none'); 
     
    end
    
    title(sprintf('%s',pt(whichPts(j)).name),'fontsize',15);
    if mod(j,10) == 1
        ylabel('Proportion of sequences');
    end
    set(gca,'fontsize',15);
  
end

for j = length(whichPts) + 1:length(pos)
    axes(ha(j));
    set(gca,'visible','off')
end

%pause
print(gcf,[destFolder,'clustBar',sprintf('%d_hours',intericTime)],'-depsc');
eps2pdf([destFolder,'clustBar',sprintf('%d_hours',intericTime),'.eps'])


%% Get clinical stuff
allOutcome = [];
allLoc = {};
for whichPt = whichPts
    [outcome(whichPt)] = getOutcome(pt,whichPt);
    allOutcome = [allOutcome;outcome(whichPt)];
    allLoc = [allLoc;pt(whichPt).clinical.seizureOnset];
    if isempty(pt(whichPt).sz_onset) == 1
        allLoc = [allLoc;nan];
    end
end

loc_bin = zeros(length(allLoc),1);
for i = 1:size(loc_bin,1)
    if isnan(allLoc{i}) == 1
        loc_bin(i) = nan;
    elseif strcmp(allLoc{i}(end-1:end),'TL') == 1
        loc_bin(i) = 1;
    else
        loc_bin(i) = 0;
    end
end

allAEDs = cell(max(whichPts),1);
onLTG = [];
for whichPt = whichPts
    allAEDs{whichPt} = getAEDs(pt(whichPt).name);
    testStr = allAEDs{whichPt};
    if isempty(testStr) == 1
        onLTG = [onLTG;nan];
    elseif any(strcmp(testStr,'LTG')) == 1
        onLTG = [onLTG;1];
    else
        onLTG = [onLTG;0];
    end
end

%% Now correlate change in pre-ic vs inter-ic with clinical stuff
% threshold p value, Bonferroni corrected
p_thresh = 0.05/length(all_p); 
differ_preic_interic = all_p < p_thresh;

[p1,info1] = correlateClinically(differ_preic_interic,allOutcome,'bin','num',0);

[p2,info2] = correlateClinically(differ_preic_interic(~isnan(loc_bin)),...
    loc_bin(~isnan(loc_bin)),'bin','bin',0);
measure = differ_preic_interic(~isnan(loc_bin));
clinical = loc_bin(~isnan(loc_bin));
figure
bar([sum(measure==0 & clinical==0) sum(measure==0 & clinical==1);...
        sum(measure==1 & clinical==0) sum(measure==1 & clinical==1)]);
xticklabels({'No preictal change','Preictal change'});
legend({'Non-temporal lobe','Temporal lobe'})
ylabel('Number of patients')


[p3,info3] = correlateClinically(differ_preic_interic(~isnan(onLTG)),...
    onLTG(~isnan(onLTG)),'bin','bin',0);

%% Now correlate change across time with clinical stuff
%{
all_p_time = [];
for whichPt = whichPts
    all_p_time = [all_p_time;p_change_time(whichPt)];
end

% threshold p value, Bonferroni corrected
p_thresh_time = 0.05/length(all_p_time); 
differ_time = all_p_time < p_thresh_time;

[p4,info4] = correlateClinically(differ_time,allOutcome,'bin','num',0);

[p5,info5] = correlateClinically(differ_time(~isnan(loc_bin)),...
    loc_bin(~isnan(loc_bin)),'bin','bin',0);
measure = differ_time(~isnan(loc_bin));
clinical = loc_bin(~isnan(loc_bin));
figure
bar([sum(measure==0 & clinical==0) sum(measure==0 & clinical==1);...
        sum(measure==1 & clinical==0) sum(measure==1 & clinical==1)]);
xticklabels({'No change across time','Change across time'});
legend({'Non-temporal lobe','Temporal lobe'})
ylabel('Number of patients')


[p6,info6] = correlateClinically(differ_time(~isnan(onLTG)),...
    onLTG(~isnan(onLTG)),'bin','bin',1);
%}

save([destFolder,'stats.mat','stats']);


end