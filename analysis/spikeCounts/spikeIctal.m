function spikeIctal(pt,cluster,whichPts)


%% Parameters

% The post-ictal time period (how many hours after the seizure I am
% defining to be post-ictal) (4 for most analyses, but 1 as a sensitivity
% analysis)
intericTime = 4;

% Plot the time periods to see what I am defining to be pre-, post-, and
% interictal. Also do some quality checks for the coverage of the
% permutation test
plotQI = 0;

plotTime = 1;

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


%% Loop through patients
for whichPt = whichPts
    
    fprintf('Doing %s\n',pt(whichPt).name);
    
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
    
    %% Get counts per bin
    counts_bin = zeros(1,size(pt(whichPt).runTimes,1));
    times_bin = zeros(1,size(pt(whichPt).runTimes,1));
    for i = 1:size(pt(whichPt).runTimes,1)
        currTimes = pt(whichPt).runTimes(i,:);
        times_bin(i) = (currTimes(1) + currTimes(2))/2;
        counts_bin(i) = sum(all_times_all >= currTimes(1) & all_times_all <= currTimes(2));
    end
    counts_all{whichPt} = counts_bin;
    times_all{whichPt} = times_bin;
    
    %% Plot times
    if plotTime == 1
        figure
        plot(times_bin/3600,counts_bin);
        hold on
        for j = 1:size(szTimes,1)
            plot([szTimes(j,1) szTimes(j,1)]/3600,ylim,'k--');
        end
        pause
        close(gcf)
    end
    
    %% Define pre-ictal, post-ictal, and interictal times
    % Define important ranges
    preIcRange = [-60*60,-1*60]; % Between 1 hour to 1 minute before a seizure
    postIcTime = 60*60*intericTime; % 60 minutes after a seizure versus 4 hours
    % Interictal range will be anything else

    preIcTimesQI = [];
    interIcTimesQI = [];
    postIcTimesQI = [];

  
    % Loop through seizures
    for j = 1:size(szTimes,1)

        % Get range of post-ictal times
        postIcTimes = szTimes(j,2) + [0 postIcTime];
        postIcTimesQI = [postIcTimesQI;postIcTimes];

      
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

    end

   

    % Get times before the first seizure
    interIcTimes(1) = min(all_times_all) + postIcTime; % assume seizure right before the record started
    interIcTimes(2) = szTimes(1,1) + preIcRange(1) - 1; % Up to 60 minutes before first sz


    if interIcTimes(1) < interIcTimes(2)
        
        interIcTimesQI = [interIcTimesQI; interIcTimes(1) interIcTimes(2)];
   
    end

    % Loop through the seizures
    for j = 1:size(szTimes,1)-1
        interIcTimes(1) = szTimes(j,2) + postIcTime;
        interIcTimes(2) = szTimes(j+1,1) + preIcRange(1) - 1;

        if interIcTimes(1) >= interIcTimes(2), continue; end
       
        interIcTimesQI = [interIcTimesQI; interIcTimes(1) interIcTimes(2)];
       
    end

    % Get time after the last seizure
    interIcTimes(1) = szTimes(end,2) + postIcTime;
    interIcTimes(2) = max(all_times_all) + preIcRange(1) -1; % Assume seizure right after record ends

    if interIcTimes(1) < interIcTimes(2)
       
        interIcTimesQI = [interIcTimesQI; interIcTimes(1) interIcTimes(2)];
        
    end
    
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
        for j = 1:size(postIcTimesQI,1)
            area([postIcTimesQI(j,1) postIcTimesQI(j,2)]/3600,[1 1],'FaceColor','b');
        end
        yl = ylim;
        for j = 1:size(szTimes,1)
           plot([szTimes(j,1) szTimes(j,1)]/3600,yl,'k--','LineWidth',5);
        end
        pause
        close(gcf)
    end
    
    %% Remove interictal times that are outside allowable
    allowableTimes = pt(whichPt).allTimes;
    interIcTimesQI = makeIntersectingTimes(interIcTimesQI,allowableTimes);
    
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
        for j = 1:size(postIcTimesQI,1)
            area([postIcTimesQI(j,1) postIcTimesQI(j,2)]/3600,[1 1],'FaceColor','b');
        end
        yl = ylim;
        for j = 1:size(szTimes,1)
           plot([szTimes(j,1) szTimes(j,1)]/3600,yl,'k--','LineWidth',5);
        end
        pause
        close(gcf)
    end
    
    %% Get numbers of spikes in appropriate times
    if 1 == 1
    num_preIc = sum(any(all_times_all' >= preIcTimesQI(:,1) & ...
        all_times_all' <= preIcTimesQI(:,2),1));
    num_interIc = sum(any(all_times_all' >= interIcTimesQI(:,1) & ...
        all_times_all' <= interIcTimesQI(:,2),1));
    num_postIc = sum(any(all_times_all' >= postIcTimesQI(:,1) & ...
        all_times_all' <= postIcTimesQI(:,2),1));
    else
    
    times_soz = all_times_all(ismember(all_spikes,soz));
    
    num_preIc = sum(any(times_soz' >= preIcTimesQI(:,1) & ...
        times_soz' <= preIcTimesQI(:,2),1));
    num_interIc = sum(any(times_soz' >= interIcTimesQI(:,1) & ...
        times_soz' <= interIcTimesQI(:,2),1));
    num_postIc = sum(any(times_soz' >= postIcTimesQI(:,1) & ...
        times_soz' <= postIcTimesQI(:,2),1));
    end
    
    %% Get spike rates
    pre_rate(whichPt) = num_preIc/sum(diff(preIcTimesQI,1,2));
    inter_rate(whichPt) = num_interIc/sum(diff(interIcTimesQI,1,2));
    post_rate(whichPt) = num_postIc/sum(diff(postIcTimesQI,1,2));
    
    %% Permutation test
    % Shuffle total number of spikes into any time period, assuming
    % randomly distributed
    %{
    p_post(whichPt) = testRates(num_interIc,sum(diff(interIcTimesQI,1,2)),...
        num_postIc,sum(diff(postIcTimesQI,1,2)),1e3);
    %}
    

    
end

pre_all = [];
inter_all = [];
post_all = [];
p_post_all = [];
temporal = [];


for whichPt = whichPts
    pre_all = [pre_all;pre_rate(whichPt)];
    inter_all = [inter_all;inter_rate(whichPt)];
    post_all = [post_all;post_rate(whichPt)];
    % p_post_all = [p_post_all;p_post(whichPt)];
    if contains(pt(whichPt).clinical.seizureOnset,'TL') == 1
        temporal = [temporal;1];
    else
        temporal = [temporal;0];
    end
end

temporal = logical(temporal);

[p,tbl,stats] = kruskalwallis([pre_all,...
inter_all,post_all],[],'off')

[p,tbl,stats] = kruskalwallis([pre_all(temporal),...
inter_all(temporal),post_all(temporal)],[],'off')


end