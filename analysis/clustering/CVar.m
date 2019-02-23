function CVar(pt,cluster,whichPts)

%{
The goal is to determine the amount of time needed to capture x% of the
total variability in cluster distribution


Plan:
- Look at most popular cluster and get the proportion in that cluster in
each hour. Then take the range from min to max proportions across all hour
long bins.
- then take all possible subsets of N consecutive hours, where I increase N such that the
width of the new range is >80% of the old range
- get range from min to max in each subset


%}

% Parameters
percVarTrue = 95; % Define true variability as nth percentiles rather than min and max
chunkMethod = 1;
doPlots = 0;
alpha = 0.8; % Try to capture what percent of variability
nboot = 1e3;

% Save file location
[~,~,~,resultsFolder,~] = fileLocations;
destFolder = [resultsFolder,'clustering/var/'];
mkdir(destFolder);

if isempty(whichPts) == 1
    for i = 1:length(pt)
        if isempty(pt(i).seq_matrix) == 0
            if size(cluster(i).bad_cluster) < cluster(i).k - 1
                whichPts = [whichPts,i];
            end
        end
    end
    if isequal(whichPts,[1,4,6,8,9,12,15,17,18,19,20,22,24,25,27,30,31]) == 0
        error('look\n');
    end
elseif whichPts == 100
    whichPts = [4 6 8 9 15 17 18 19 20 22 24 25 27 30 31];
end

allMinCapture = [];
allMinProp = [];
allHours = [];


for whichPt = whichPts
    
    fprintf('Doing %s\n',pt(whichPt).name);
    
    % Get patient parameters
    locs = pt(whichPt).electrodeData.locs(:,2:4);
    szTimes = pt(whichPt).newSzTimes;
    soz = pt(whichPt).newSOZChs;
    saveFolder = [destFolder,pt(whichPt).name,'/'];
    mkdir(saveFolder);
    
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
    
    % Get most popular cluster
    
    
    %% Get actual time chunks
    test_t = 3600; % 60 minute chunks
    allTimes = pt(whichPt).allTimes;
    
    if chunkMethod == 1
        % New way
        [prop_pop,chunk_times,~] = getChunkTimes(allTimes,test_t,all_times_all,idx);
        n_chunks = length(chunk_times);
    else
    % Old way
        popular = mode(idx);
    
        % Divide run into 60 minute chunks
        % This may result in some chunks that are empty because there was low
        % voltage data and so I skipped spike detection in this period, but
        % this should not affect the cluster distribution. 
        n_chunks = ceil((max(all_times_all) - min(all_times_all))/test_t);
        prop_pop = zeros(n_chunks,1);
        chunk_times = zeros(n_chunks,1);

        for i = 1:n_chunks

            % Get the time range for the chunk
            curr_times = [min(all_times_all) + (i-1)*test_t,...
               min(min(all_times_all) + i*test_t,max(all_times_all))];

            % Get the spike indices in that time chunk
            chunk_spikes = find(all_times_all >= curr_times(1) & ...
                all_times_all <= curr_times(2));

            chunk_times(i) = curr_times(1);
            
            % Get the proportion of spikes in the most popular cluster
            prop_pop(i) = sum(idx(chunk_spikes) == popular)/length(chunk_spikes);

        end
    end

    %{
    figure
    plot(chunk_times,prop_pop)
    %}
    
    
    true_range_old = [min(prop_pop) max(prop_pop)]
    true_range = [prctile(prop_pop,100-percVarTrue),prctile(prop_pop,percVarTrue)]
    
    %% Now, start with 1 hour and go up, and get the new range
    sub_range = zeros(n_chunks,2);
    sub_std = zeros(n_chunks,2);
    for i = 1:n_chunks
        
        
        
        % Get all subsets of n continuous hours
        allSubsets = getContinuousSubset(n_chunks,i);
        
        % Need to re-initialize temp_range
        temp_range = zeros(size(allSubsets,1),2);
        
        for j = 1:size(allSubsets,1)
            subset = allSubsets(j,:);
            temp_prop_pop = prop_pop(subset);
            temp_range(j,:) = [prctile(temp_prop_pop,100-percVarTrue),prctile(temp_prop_pop,percVarTrue)];
        end
        sub_range(i,:) = mean(temp_range,1);
        
        if i == n_chunks
            sub_std(i,:) = [0 0];
        else
            sub_std(i,:) = std(temp_range,1);
        end
        
        
        
        % Old way: 1000 random non-continuous subsets
        %fprintf('Doing %d of %d\n',i,n_chunks)
        %{
        % Do 1000 trials where I randomly select one of those combos
        temp_range = zeros(nboot,2);
        for ib = 1:nboot
            subset = randsample(n_chunks,i);
            temp_prop_pop = prop_pop(subset);
            temp_range(ib,:) = [min(temp_prop_pop) max(temp_prop_pop)];
        end
        sub_range(i,:) = mean(temp_range,1);
        sub_std(i,:) = std(temp_range,1);
        %}
        
    end
    
    true_width = diff(true_range);
    sub_width = diff(sub_range,1,2);
    
    prop_true_width = sub_width/true_width;
    capture_var = find(prop_true_width>alpha);
    min_capture_var = min(capture_var);
    
    allMinCapture = [allMinCapture; min_capture_var];
    allMinProp = [allMinProp; min_capture_var/n_chunks];
    allHours = [allHours;n_chunks];
    
    outcome(whichPt) = getOutcome(pt,whichPt);
    
    %% Plot
    if doPlots == 1
        figure
        set(gcf,'Position',[440 233 977 565])
        mn = errorbar(1:n_chunks,sub_range(:,1),sub_std(:,1),'k','LineWidth',2);
        hold on
        mx = errorbar(1:n_chunks,sub_range(:,2),sub_std(:,2),'k','LineWidth',2);


        xPos = [min_capture_var+10 min_capture_var];
        yPos = [mean(true_range), mean(true_range)];

        pos = get(gca,'Position');
        annotation('textarrow',[(xPos(1) + abs(min(xlim)))/diff(xlim) * pos(3) + pos(1),...
            (xPos(2) + abs(min(xlim)))/diff(xlim) * pos(3) + pos(1)],...
            [(yPos(1) - (min(ylim)))/diff(ylim) * pos(4) + pos(2),...
            (yPos(2) - (min(ylim)))/diff(ylim) * pos(4) + pos(2)],...
            'String',...
            sprintf('Duration needed to capture\n%d%% of variability: %d hours\n(%1.1f%% of total dataset)', ...
            alpha*100,min_capture_var,min_capture_var/n_chunks*100),'FontSize',25);

        xlabel('Number of consecutive hours');
        ylabel(sprintf('5th and 95th %%ile proportion of spikes in\nmost popular cluster across hour-long bins'));
        title(sprintf('Dependence on sampling of\nspike location variability for %s',...
            pt(whichPt).name),'fontsize',25);
        set(gca,'fontsize',25);
        hold on
        cp = plot([min_capture_var, min_capture_var],get(gca,'ylim'),'k--','LineWidth',3);
        cp = plot([min_capture_var, min_capture_var],get(gca,'ylim'),'k--','LineWidth',3); % needed for voodoo
       
        ylim([min(ylim),min(max(ylim),1.2)])
        pause
        print(gcf,[saveFolder,'clustVar_',sprintf('%d_%s',percVarTrue,pt(whichPt).name)],'-depsc');
        eps2pdf([saveFolder,'clustVar_',sprintf('%d_%s',percVarTrue,pt(whichPt).name),'.eps'])
        close(gcf)
    end
    
end

fprintf('The mean number of hours needed to capture 80%% of the variability for %d%% was %1.1f (range %d-%d)\n',...
    percVarTrue,mean(allMinCapture),min(allMinCapture),max(allMinCapture));

fprintf('The mean percent of hours needed to capture 80%% of the variability for %d%% was %1.1f (range %1.1f-%1.1f)\n',...
    percVarTrue,mean(allMinProp)*100,min(allMinProp)*100,max(allMinProp)*100);

%% Do stats correlating number needed to capture with outcome

allOutcome = [];
allLoc = {};
for whichPt = whichPts
    allOutcome = [allOutcome;outcome(whichPt)];
    allLoc = [allLoc;pt(whichPt).clinical.seizureOnset];
    %{
    if isempty(pt(whichPt).sz_onset) == 1
        allLoc = [allLoc;nan];
    end
    %}
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

[p1,info1] = correlateClinically(allMinCapture,allOutcome,'num','num',0);
fprintf(['The p-value for Spearman rank correlation of number of hours needed\n'...
    'to capture variability and outcome is p = %1.2e.\n'],p1);


[p2,info2] = correlateClinically(allMinCapture(~isnan(loc_bin)),...
loc_bin(~isnan(loc_bin)),'num','bin',0);
fprintf(['The p-value for Wilcoxon rank sum of number of hours needed\n'...
    'to capture variability and temporal vs non temporal lobe is p = %1.2e.\n'],p2);

% Get minimum hours needed to capture for temporal and extra-temporal
% seizures
min_capture_temp = allMinCapture(loc_bin ==1);
min_capture_not_temp = allMinCapture(loc_bin ~= 1);

fprintf('Mean hours for temp is %1.1f and for extra-temp is %1.1f\n',...
    median(min_capture_temp),median(min_capture_not_temp));

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
[p3,info3] = correlateClinically(allMinCapture(~isnan(onLTG)),...
    onLTG(~isnan(onLTG)),'num','bin',0);

% Spearman correlation coefficient (non parametric rank)
%{
[rho,pval] = corr(allMinCapture,allOutcome,...
    'Type','Spearman');

%% Plot
figure
scatter(allMinCapture,allOutcome,100,'filled');
xlabel('Minimum number of hours needed to capture 80 percent variance');
ylabel('Outcome (modified Engel)');
title(sprintf('Hours needed to capture variance vs outcome, p = %1.2f',pval))


[~,pval2] = corr(allHours,allOutcome,...
    'Type','Spearman');
figure
scatter(allHours,allOutcome,100,'filled');
xlabel('Total number of hours');
ylabel('Outcome (modified Engel)');
title(sprintf('Hours vs outcome, p = %1.2f',pval2))
%}

%allMinCapture
%allMinProp

end