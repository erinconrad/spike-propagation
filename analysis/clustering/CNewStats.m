function CNewStats(pt,cluster,whichPts)

%{ 

CStats
This is my cleaned up file for getting statistics on the cluster data 

%}

% Parameters
intericTime = 4;

plotQI = 1;
doPermPlot = 0;
doPlots = 1;
doLongStuff = 1;
doPre = 1;

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



allCounts = [];
allPat = [];
allChunk = [];
chi_tables_plot = cell(max(whichPts),1);
p_plot = zeros(max(whichPts),1);
p_change_time = zeros(max(whichPts),1);
allDist = [];

chi_tables_postIc = cell(max(whichPts),1);
p_postIc = zeros(max(whichPts),1);

preIcDistAll = [];
interIcDistAll = [];
postIcDistAll = [];

preIcDistMean = [];
interIcDistMean = [];
postIcDistMean = [];

diffDistAll = [];
pDistAll = [];

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
    
    if isempty(soz) == 1
        spike_dist = [];
    else
        spike_dist = zeros(size(all_locs,1),1);
        for i = 1:size(all_locs,1)
            spike_dist(i) = min(vecnorm(all_locs(i,:) - locs(soz,:),2,2)); 
        end
    end
    
    % Analysis 0
    dist = vecnorm(diff(C,1),2,2);
    allDist = [allDist;dist];
    
    
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
    
    if isnan(p_1) == 1
        if k == length(bad_cluster) + 1
            p_1 = 1;
            fprintf('Only one cluster, defining p-value to be 1\n');
        else
            error('What\n');
        end
    end
    
    % Save information into patient struct
    stats(whichPt).hour.tbl = tbl_1;
    stats(whichPt).hour.chi2 = chi2_1;
    stats(whichPt).hour.p = p_1;
    stats(whichPt).hour.labels = labels_1;
    
    fprintf(['For %s, regarding whether 60 minute chunks\n have different cluster'...
    ' distributions,\n the p-value is %1.1e\n\n\n'],pt(whichPt).name,p_1);
    
    p_change_time(whichPt) = p_1;
    
    if doLongStuff == 1
    
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
        postIcTimesQI = [];

        % Get all the pre-ictal and post-ictal spikes
        preIcSpikes = [];
        postIcSpikes = [];
        preIcSpikeTimes = [];
        postIcSpikeTimes = [];
        preIcSpikeNums = [];
        postIcSpikeNums = [];

        % Loop through seizures
        for j = 1:size(szTimes,1)

            % Get range of post-ictal times
            postIcTimes = szTimes(j,2) + [0 postIcTime];
            postIcTimesQI = [postIcTimesQI;postIcTimes];

            % Get the indices of spikes in this time range
            post_ic_spike_idx = (all_times_all >= postIcTimes(1) & ...
                all_times_all <= postIcTimes(2));
            s = find(post_ic_spike_idx);
            % Remove those that we've already included
            s(ismember(s,postIcSpikes)) = [];
            postIcSpikes = [postIcSpikes;s];
            postIcSpikeTimes = [postIcSpikeTimes;all_times_all(s)];
            postIcSpikeNums = [postIcSpikeNums;length(s)];
            
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

            % Get the cluster indices, times of spikes, and number of
            % spikes
            preIcSpikes = [preIcSpikes;find(spike_idx)];
            preIcSpikeTimes = [preIcSpikeTimes;all_times_all(spike_idx)];
            preIcSpikeNums = [preIcSpikeNums;sum(spike_idx)];

        end

        % Get all the inter-ictal spikes
        interIcSpikes = [];
        interIcSpikeTimes = [];
        interIcSpikeNums = [];

        % Get times before the first seizure
        interIcTimes(1) = min(all_times_all) + postIcTime; % assume seizure right before the record started
        interIcTimes(2) = szTimes(1,1) + preIcRange(1) - 1; % Up to 60 minutes before first sz


        if interIcTimes(1) < interIcTimes(2)
            spike_idx = all_times_all >= interIcTimes(1) & ...
                all_times_all <= interIcTimes(2);
            interIcTimesQI = [interIcTimesQI; interIcTimes(1) interIcTimes(2)];
            interIcSpikes = [interIcSpikes;find(spike_idx)];
            interIcSpikeTimes = [interIcSpikeTimes;all_times_all(spike_idx)];
            interIcSpikeNums = [interIcSpikeNums;sum(spike_idx)];
        end

        % Loop through the seizures
        for j = 1:size(szTimes,1)-1
            interIcTimes(1) = szTimes(j,2) + postIcTime;
            interIcTimes(2) = szTimes(j+1,1) + preIcRange(1) - 1;

            if interIcTimes(1) >= interIcTimes(2), continue; end
            spike_idx = all_times_all >= interIcTimes(1) & ...
                all_times_all <= interIcTimes(2);
            interIcTimesQI = [interIcTimesQI; interIcTimes(1) interIcTimes(2)];
            interIcSpikes = [interIcSpikes;find(spike_idx)];
            interIcSpikeTimes = [interIcSpikeTimes;all_times_all(spike_idx)];
            interIcSpikeNums = [interIcSpikeNums;sum(spike_idx)];
        end

        % Get time after the last seizure
        interIcTimes(1) = szTimes(end,2) + postIcTime;
        interIcTimes(2) = max(all_times_all) + preIcRange(1) -1; % Assume seizure right after record ends

        if interIcTimes(1) < interIcTimes(2)
            spike_idx = all_times_all >= interIcTimes(1) & ...
                all_times_all <= interIcTimes(2);
            interIcTimesQI = [interIcTimesQI; interIcTimes(1) interIcTimes(2)];
            interIcSpikes = [interIcSpikes;find(spike_idx)];
            interIcSpikeTimes = [interIcSpikeTimes;all_times_all(spike_idx)];
            interIcSpikeNums = [interIcSpikeNums;sum(spike_idx)];
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

        

       %{
        Analysis 2: Is the cluster distribution different in the pre-ictal
        compared to the interictal period?

       %}
        
        if doPre == 1
        all_s = [(preIcSpikes);(interIcSpikes)];
        all_s_times = [preIcSpikeTimes;interIcSpikeTimes];
        colors = [repmat([1,0,0],length(preIcSpikes),1);...
            repmat([0,0,1],length(interIcSpikes),1)];
        all_s_old = all_s;
        [B,I] = sort(all_s_times);
        all_s = all_s(I);

        % To test it, show the times of all the sorted pre and interictal
        % spikes and their cluster distributions
        if doPermPlot == 1
            new_colors = colors(I,:);
            figure
            scatter(B/3600,ones(length(B),1),50,new_colors);
            hold on
            szTimes = pt(whichPt).newSzTimes;
            for j = 1:size(szTimes,1)
                plot([szTimes(j,1),szTimes(j,1)]/3600,ylim,'k--');
            end
            pause
            close(gcf)
        end
        
        % Now, for the permutation test, I will take random chunks of
        % spikes equal in number to the original pre-ictal time periods
        nboot = 1e3;
        n_spikes = length(all_s);
        chi2_boot = zeros(nboot,1);
        for ib = 1:nboot
            
            if mod(ib,100) == 0
                fprintf('Doing %d of %d\n',ib,nboot);
            end
            
            % The number of spikes in each pre-ictal period
            new_pre = cell(size(preIcSpikeNums,1),1);
            
            % Loop through pre-ictal periods
            for j = 1:size(preIcSpikeNums,1)
                
                % How many spikes to pick (number of true pre-ictal spikes
                % in that period)
                n_chunk = preIcSpikeNums(j);
                
                % Try to generate a random set of spikes equal to that
                % number
                while 1
                    % pick a start index
                    start = randi(n_spikes);

                    % See if it will run into one of the other chunks
                    if isempty(intersect(mod(start-1:start+n_chunk-1,n_spikes) + 1,...
                            unique([new_pre{:}]))) == 1
                        break
                    end
                end
                
                % If it doesn't run into one of the other chunks, add it to
                % the new fake pre-ictal spikes
                new_pre{j} = mod(start-1:start+n_chunk-1,n_spikes) + 1;
                
            end
            
           % There should not be ANY repeats
           if isequal(sort([new_pre{:}]),unique([new_pre{:}])) == 0
               error('What\n');
           end
            
            % Get interictal indices
            new_inter = setdiff(1:length(all_s),unique([new_pre{:}]));
            
            % Get the actual spike indices
            pre_s = all_s(unique([new_pre{:}]));
            inter_s = all_s(new_inter);
            
            % Get their clusters
            pre_clust = idx(pre_s);
            inter_clust = idx(inter_s);
            
            if 1==0
                % test this by plotting the times of these fake pre-ictal
                % and interictal spikes
                figure
                scatter(all_times_all(inter_s)/3600,ones(length(inter_s),1),50,'r');
                hold on
                scatter(all_times_all(pre_s)/3600,ones(length(pre_s),1),50,'b');


                for j = 1:size(szTimes,1)
                    plot([szTimes(j,1),szTimes(j,1)]/3600,ylim,'k--');
                end
                text(mean(xlim),[1.5],sprintf('%d pre-ic spikes',length(pre_s)));
                pause
                close(gcf)
            end
            
            [~,chi2_boot(ib)] = crosstab([ones(length(pre_clust),1);...
                2*ones(length(inter_clust),1)],[pre_clust;inter_clust]);
            
        end

       

        % Do it once for real
        inter_clust = idx(interIcSpikes);
        pre_clust = idx(preIcSpikes);
        [tbl_2,chi2_real] = crosstab([ones(length(pre_clust),1);...
                2*ones(length(inter_clust),1)],[pre_clust;inter_clust]);

        % Sort the permutation test chi squareds and find how many are
        % larger than the real chi squared
        sorted_boot = sort(chi2_boot);
        diff_s_chi = sorted_boot-chi2_real;
        allLarger = find(diff_s_chi>0);
        
        % The p value is the percentage larger 
        if isempty(allLarger) == 1
            p_2 = 0;
        else
            firstLarger = allLarger(1);
            p_2 = (nboot-firstLarger+1)/(nboot+1);
        end

        % Plot the result of the permuation test
        if doPermPlot == 1
            figure 
            scatter(1:length(sorted_boot),sorted_boot);
            hold on
            plot(xlim,[chi2_real chi2_real]);
            text(mean(xlim),mean(ylim),sprintf('%1.1e',p_2));

            pause
            close(gcf)
        end
        %}

        if k == length(bad_cluster) + 1
            p_2 = 1;
            fprintf('Only one cluster, defining p-value to be 1\n');
        end


        % Save information 
        stats(whichPt).preIc.tbl = tbl_2;
        stats(whichPt).preIc.p = p_2;
        chi_tables_plot{whichPt} = tbl_2;
        p_plot(whichPt) = p_2;
        end


        %{
        Analysis 3: post-ictal vs the rest
        %}
        all_s = [(preIcSpikes);(interIcSpikes);postIcSpikes];
        all_s_times = [preIcSpikeTimes;interIcSpikeTimes;postIcSpikeTimes];
        colors = [repmat([1,0,0],length(preIcSpikes),1);...
            repmat([0,0,1],length(interIcSpikes),1);...
            repmat([0,1,0],length(postIcSpikes),1)];

        [B,I] = sort(all_s_times);
        all_s = all_s(I);

        if doPermPlot == 1
            new_colors = colors(I,:);
            figure
            scatter(B/3600,all_s,50,new_colors);
            hold on
            szTimes = pt(whichPt).newSzTimes;
            for j = 1:size(szTimes,1)
                plot([szTimes(j,1),szTimes(j,1)]/3600,ylim,'k--');
            end
            pause
            close(gcf)
        end
        
        if length(postIcSpikes) ~= length(unique(postIcSpikes))
            error('What\n');
        end
        
        nboot = 1e3;
        n_spikes = length(all_s);
        chi2_boot = zeros(nboot,1);
        dist_diff_boot = zeros(nboot,1);
        
        % Sort the post ic spike numbers in descending order. I do this
        % because I am trying to fit a bunch of smaller chunks, non
        % overlapping, into a big chunk. If I put in the small chunks
        % first, I am likely to put them in a configuration that will not
        % allow me to fit the bigger ones in later, so my while loop gets
        % stuck!!!
        postIcSpikeNums = sort(postIcSpikeNums,'descend');
        for ib = 1:nboot
            
            if mod(ib,10) == 0
                fprintf('Doing %d of %d\n',ib,nboot);
            end
            
                new_post = cell(size(postIcSpikeNums,1),1);

                %Get post-ictal indices
                for j = 1:size(postIcSpikeNums,1)
                    
                    % How many spikes to pick
                    n_chunk = postIcSpikeNums(j);


                    while 1
                        % pick a start index
                        start = randi(n_spikes);

                        % See if it will run into one of the other chunks
                        if isempty(intersect(mod(start-1:start+n_chunk-1,n_spikes) + 1,...
                                unique([new_post{:}]))) == 1
                            break
                        end
                        
                       
                    end
                    
                    new_post{j} = mod(start-1:start+n_chunk-1,n_spikes) + 1;                 
                end
                
                
            
           if isequal(sort([new_post{:}]),unique([new_post{:}])) == 0
               error('What\n');
           end
           
           if sum(postIcSpikeNums) ~= length(postIcSpikes)
               error('what\n');
           end
            
            % Get other indices
            new_other = setdiff(1:length(all_s),unique([new_post{:}]));
            post_s = all_s(unique([new_post{:}]));
            other_s = all_s(new_other);
            post_clust = idx(post_s);
            other_clust = idx(other_s);
            
            % Get distances
            if isempty(soz) == 1
                diff_dist_boot(ib) = nan;
            else
                boot_post_dist = mean(spike_dist(post_s));
                boot_other_dist = mean(spike_dist(other_s));
                diff_dist_boot(ib) = boot_post_dist-boot_other_dist;
            end
            
            if 1== 0
            figure
            scatter(all_times_all(other_s)/3600,ones(length(other_s),1),50,'r');
            hold on
            scatter(all_times_all(post_s)/3600,ones(length(post_s),1),50,'g');
  
            
            for j = 1:size(szTimes,1)
                plot([szTimes(j,1),szTimes(j,1)]/3600,ylim,'k--');
            end
            text(mean(xlim),[1.5],sprintf('%d post-ic spikes',length(post_s)));
            pause
            close(gcf)
            end
            
            [~,chi2_boot(ib)] = crosstab([ones(length(post_clust),1);...
                2*ones(length(other_clust),1)],[post_clust;other_clust]);
            
        end
  
        % Do it once for real
        post_clust = idx(postIcSpikes);
        other_clust = idx([preIcSpikes;interIcSpikes]);
        [tbl_3,chi2_real] = crosstab([ones(length(post_clust),1);...
                2*ones(length(other_clust),1)],[post_clust;other_clust]);


        sorted_boot = sort(chi2_boot);

        diff_s_chi = sorted_boot-chi2_real;
        allLarger = find(diff_s_chi>0);
        if isempty(allLarger) == 1
            p_3 = 0;
        else
            firstLarger = allLarger(1);
            p_3 = (nboot-firstLarger)/nboot;
        end

        if doPermPlot == 1
            figure 
            scatter(1:length(sorted_boot),sorted_boot);
            hold on
            plot(xlim,[chi2_real chi2_real]);
            text(mean(xlim),mean(ylim),sprintf('%1.1e',p_3));

            pause
            close(gcf)
        end
        %}



        if k == length(bad_cluster) + 1
            p_3 = 1;
            fprintf('Only one cluster, defining p-value to be 1\n');
        end


        % Save information into patient struct
        stats(whichPt).postIc.tbl = tbl_3;
        stats(whichPt).postIc.p = p_3;
        chi_tables_postIc{whichPt} = tbl_3;
        p_postIc(whichPt) = p_3;
        
        
        %% Distance analysis
        % Do it once for real
        if isempty(soz) == 1
            p_dist = [];
            sorted_boot = [];
            diff_dist_real = [];
        else
            post_dist = mean(spike_dist(postIcSpikes));
            other_dist = mean(spike_dist([preIcSpikes;interIcSpikes]));
            diff_dist_real = post_dist-other_dist;

            sorted_boot = sort(diff_dist_boot);
            p_dist = (sum(abs(sorted_boot) > abs(diff_dist_real))+1)/...
                (length(sorted_boot) + 1);
        end
        
        pDistAll = [pDistAll;p_dist];
        diffDistAll = [diffDistAll;diff_dist_real];
        
        stats(whichPt).dist.p = p_dist;
        stats(whichPt).dist.dist_diff = diff_dist_real;
        
        if doPermPlot == 1 && isempty(soz) == 0
            figure
            scatter(1:length(sorted_boot),sorted_boot)
            hold on
            plot(xlim,[diff_dist_real diff_dist_real]);
            text(mean(xlim),mean(ylim),sprintf('%1.1e',p_dist));
            pause
            close(gcf)
            
        end
        
        
        
    end
    
   
    
end

if intericTime == 4
    save([destFolder,'stats4.mat'],'stats');
elseif intericTime == 1
    save([destFolder,'stats1.mat'],'stats');
end

%% Fisher's method to combine p values for change over time
all_p_change = [];
for whichPt = whichPts
    all_p_change= [all_p_change;p_change_time(whichPt)];
end
X_2 = -2 * sum(log(all_p_change));
sum_p = 1-chi2cdf(X_2,2*length(all_p_change));
fprintf('%d out of %d patients showed a significant change in cluster distribution over time.\n',...
    sum(all_p_change < 0.05/length(whichPts)),length(whichPts))
fprintf('The group p value for change over time is %1.1e\n',sum_p);

if doLongStuff == 1
    
    if doPre == 1
    %% Fisher's method to combine p values for pre-ic vs interic
    all_p = [];
    for whichPt = whichPts
       all_p = [all_p;p_plot(whichPt)]; 
    end

    X_2 = -2 * sum(log(all_p));
    sum_p = 1-chi2cdf(X_2,2*length(all_p));

    fprintf('The group p value for pre-ictal vs interictal cluster is %1.1e\n',sum_p);
    fprintf('%d of %d patients had a different preictal vs interictal distribution.\n',...
        sum(all_p<0.05/length(whichPts)),length(whichPts));
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
    end

    %% Bar graphs
    if doPlots == 1
        if doPre == 1
        %% Pre-ic
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
            xticklabels({'Pre','Inter'})

            % Plot p-value
            p = p_plot(whichPts(j));
            if p < 0.001/length(whichPts)
                flag = '***';
            elseif p < 0.01/length(whichPts)
                flag = '**';
            elseif p < 0.05/length(whichPts)
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

            if j == 17
               lgnd=legend(b,legend_names,'Position',...
                   [pos{j}(1)+0.025 pos{j}(2)+0.25,...
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
        end

        %% Post-ic
        
        all_p = [];
        for whichPt = whichPts
           all_p = [all_p;p_postIc(whichPt)]; 
        end

        X_2 = -2 * sum(log(all_p));
        sum_p = 1-chi2cdf(X_2,2*length(all_p));

        fprintf('The group p value for pre-ictal vs interictal cluster is %1.1e\n',sum_p);
        fprintf('%d of %d patients had a different preictal vs interictal distribution.\n',...
            sum(all_p<0.05/length(whichPts)),length(whichPts));
        
        figure
        set(gcf,'Position',[99 26 1264 686]);
        numcols = 10;
        numrows = ceil(length(whichPts)/numcols);
        [ha, pos] = tight_subplot(numrows,numcols,[.08 .01],[.05 .05],[.05 .01]); 
        for j = 1:length(whichPts)
            axes(ha(j));

            % Get chi2 table
            tbl = chi_tables_postIc{whichPts(j)};

            % Get proportions rather than absolute numbers
            prop = tbl./sum(tbl,2);

            b=bar(prop);
            xticklabels({'Post','Other'})

            % Plot p-value
            p = p_postIc(whichPts(j));
            if p < 0.001/length(whichPts)
                flag = '***';
            elseif p < 0.01/length(whichPts)
                flag = '**';
            elseif p < 0.05/length(whichPts)
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

            if j == 17
               lgnd=legend(b,legend_names,'Position',...
                   [pos{j}(1)+0.025 pos{j}(2)+0.25,...
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
        print(gcf,[destFolder,'clustPost',sprintf('%d_hours',intericTime)],'-depsc');
        eps2pdf([destFolder,'clustPost',sprintf('%d_hours',intericTime),'.eps'])


    end
   

end


end