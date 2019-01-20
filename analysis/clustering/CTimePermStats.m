function CTimePermStats(pt,cluster,whichPts)

%{ 

CStats
This is my cleaned up file for getting statistics on the cluster data 

%}

% Parameters
intericTime = 4;

plotQI = 0;
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


        % Loop through seizures
        for j = 1:size(szTimes,1)

            % Get range of post-ictal times
            postIcTimes = szTimes(j,2) + [0 postIcTime];
            for kk = 1:size(postIcTimesQI,1)
                % shorten it if in range from prior postIctalTime
                if postIcTimesQI(kk,2) > postIcTimes(1)
                    postIcTimesQI(kk,2) = postIcTimes(2);
                    postIcTimes = [];
                    break
                end
            end
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
            spike_idx = all_times_all >= interIcTimes(1) & ...
                all_times_all <= interIcTimes(2);
            interIcTimesQI = [interIcTimesQI; interIcTimes(1) interIcTimes(2)];
        end

        % Loop through the seizures
        for j = 1:size(szTimes,1)-1
            interIcTimes(1) = szTimes(j,2) + postIcTime;
            interIcTimes(2) = szTimes(j+1,1) + preIcRange(1) - 1;

            if interIcTimes(1) >= interIcTimes(2), continue; end
            spike_idx = all_times_all >= interIcTimes(1) & ...
                all_times_all <= interIcTimes(2);
            interIcTimesQI = [interIcTimesQI; interIcTimes(1) interIcTimes(2)];
           
        end

        % Get time after the last seizure
        interIcTimes(1) = szTimes(end,2) + postIcTime;
        interIcTimes(2) = max(all_times_all) + preIcRange(1) -1; % Assume seizure right after record ends

        if interIcTimes(1) < interIcTimes(2)
            spike_idx = all_times_all >= interIcTimes(1) & ...
                all_times_all <= interIcTimes(2);
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

        
        %% Pre-ic analysis
        
        % Concatenate all continuous blocks of time. Very commonly a
        % preictal block will be immediately be preceded by an interictal
        % block. I want to smush those together
        allTimes = [interIcTimesQI;preIcTimesQI];
        allTimes = sortrows(allTimes,1);
        allowableBlocks = allTimes(1,:);
        for i = 2:size(allTimes,1)
            % If the first time here is very close to the last time of the
            % last allowable block, concatenate them
            if abs(allTimes(i,1) - allowableBlocks(end,2)) < 60 % 60 seconds apart
                allowableBlocks(end,2) = allTimes(i,2);
                
            % If not, add this row as a next chunk of allowable time
            else
                allowableBlocks = [allowableBlocks;allTimes(i,:)];
            end
        end
        
        % I now have an array of all allowable chunks of time. My plan will
        % be for each permutation to pick an integer going from 1 to
        % allTotalTime. I then find which cumsum it is in and that tells me
        % which block to try to stick it in.
        totalTimeInChunks = diff(allowableBlocks,1,2);
        allTotalTime = round(sum(totalTimeInChunks));
        cumSumTimes = cumsum(totalTimeInChunks);
        
        
        % Sort preIcTimes by duration (doesn't matter for preIcTimes but
        % will when I replicate this for postIcTimes).
        preIcTimesDuration = diff(preIcTimesQI,1,2);
        [preIcTimesDurationSorted,I] = sort(preIcTimesDuration,'descend');
        preIcTimesSorted = preIcTimesQI(I,:);
        
        nboot = 1e3;
        chi2_boot = zeros(nboot,1);
        
        for ib = 1:nboot
            
            if mod(ib,100) == 0
                fprintf('Doing %d of %d\n',ib,nboot);
            end
            
            fakePreIcTimes = zeros(size(preIcTimesSorted));
            
            % which pre-ictal time am I faking?
            for j = 1:size(fakePreIcTimes,1)
                
                while 1
            
                    t = randi(allTotalTime);

                    % Find which cumsum t fits into. It fits into the last
                    % non-negative one
                    difft = find(t-cumSumTimes>0);
                    if isempty(difft) == 1
                        whichBlock = 1;
                    else
                        whichBlock = difft(end)+1;
                    end
                    allowableTimes = allowableBlocks(whichBlock,:);
                    blockDur = diff(allowableTimes);

                    % Now I know which allowable block to stick it in. Pick
                    % a random time in that block, far enough from the end.
                    currDur = preIcTimesDurationSorted(j); % Duration of current preictal block
                    timeLeft = blockDur -  currDur;
                    if timeLeft < 0
                        continue
                    end
                    t2 = randi(round(timeLeft+1));
                    
                    % The candidate fake pre-ictal time
                    fakePreIcTime = [allowableTimes(1) + t2,...
                        allowableTimes(1) + t2 + currDur];
                    
                    % Throw it out and try again if for some reason it goes
                    % beyond the end of the block (it shouldn't)
                    if fakePreIcTime(2) > allowableTimes(2)
                        continue
                    end
                    
                    intersect = 0;
                    
                    % Throw it out and try again if it intersects any of
                    % the previous fake pre-ictal times
                    for kk = 1:j-1
                        oldRange = fakePreIcTimes(kk,:);
                        
                        % If intersect
                        if doTimeRangesIntersect(fakePreIcTime,oldRange) == 1
                            intersect = 1;
                            break
                        end
                        
                    end
                    
                    if intersect == 1, continue; end
                    
                    % If it makes to here, it is acceptable
                    fakePreIcTimes(j,:) = fakePreIcTime;
                    break;
                
                end
                
            end
            
            % Resort the fake preIctal times
            fakePreIcTimes = sortrows(fakePreIcTimes);
            
            % Do a check for intersections
            intersection = 0;
            
            for i = 1:size(fakePreIcTimes)
                for j = i+1:size(fakePreIcTimes)
                    if doTimeRangesIntersect(fakePreIcTimes(i,:),...
                            fakePreIcTimes(j,:)) == 1
                        intersection = 1;
                    end
                end
            end
            
            if intersection == 1
                error('What\n');
            end
           
            % Now I need to build the fake interictal times
            fakeInterIcTimes = [];
            
            % Start with the first block and the first pre-ictal time
            startTime = allowableBlocks(1,1);
            whichBlock = 1;
            whichPreIc = 1;
            
            while 1
                currRange = [startTime allowableBlocks(whichBlock,2)];
                if doTimeRangesIntersect(currRange,fakePreIcTimes(whichPreIc,:)) == 1
                    fakeInterIcTimes = [fakeInterIcTimes;...
                        startTime fakePreIcTimes(whichPreIc,1)];
                    startTime = fakePreIcTimes(whichPreIc,2);
                    whichPreIc = whichPreIc + 1;
                else
                    whichPreIc = whichPreIc + 1;
                end
                
                % If I've reached the end of the pre-ictal blocks
                if whichPreIc == size(fakePreIcTimes,1) + 1
                    
                    % Add the remaining times
                    fakeInterIcTimes = [fakeInterIcTimes;...
                            startTime allowableBlocks(whichBlock,2)];
                    
                    % Increase the block and start again
                    whichBlock = whichBlock + 1;
                    if whichBlock == size(allowableBlocks,1) + 1
                        break
                    end
                    startTime = allowableBlocks(whichBlock,1);
                    whichPreIc = 1;
                        
                end 
            end
            
            % Do a check for intersections
            intersection = 0;
            
            for i = 1:size(fakeInterIcTimes)
                for j = i+1:size(fakeInterIcTimes)
                    if doTimeRangesIntersect(fakeInterIcTimes(i,:),...
                            fakeInterIcTimes(j,:)) == 1
                        intersection = 1;
                    end
                end
            end
            
            if intersection == 1
                error('What\n');
            end
            
            % Do a final check for any intersections
            intersection = 0;
            
            for i = 1:size(fakeInterIcTimes)
                for j = 1:size(fakePreIcTimes)
                    if doTimeRangesIntersectForgiving(fakeInterIcTimes(i,:),...
                            fakePreIcTimes(j,:)) == 1
                        intersection = 1;
                    end
                end
            end
            
            if intersection == 1
                error('What\n');
            end
            
            % Do a final check to make sure times add up
            totalTimeFake = sum(diff(fakeInterIcTimes,1,2)) + ...
                sum(diff(fakePreIcTimes,1,2));
            if abs(totalTimeFake-allTotalTime) > 2
                error('What\n');
            end
            
            
            % Plot the new times
            if 1 == 0
                
                f=figure;
                
                for j = 1:size(postIcTimesQI,1)
                    area([postIcTimesQI(j,1) postIcTimesQI(j,2)]/3600,[1 1],'FaceColor','b');
                    hold on
                end
                
                
                
                for j = 1:size(fakePreIcTimes,1)
                   area([fakePreIcTimes(j,1) fakePreIcTimes(j,2)]/3600,[1 1],'FaceColor','g');
                   hold on
                end
                
                
                for j = 1:size(fakeInterIcTimes,1)
                   area([fakeInterIcTimes(j,1) fakeInterIcTimes(j,2)]/3600,[1 1],'FaceColor','r');
                   hold on
                end
                
                yl = ylim;
                for j = 1:size(szTimes,1)
                   plot([szTimes(j,1) szTimes(j,1)]/3600,yl,'k--','LineWidth',5);
                end
                pause
                close(f);
                
            end
            
            % Now I have the correct times, let me get the spikes in the
            % times
            fakePreIcSpikes = [];
            for i = 1:size(fakePreIcTimes,1)
                spike_idx = find(all_times_all >= fakePreIcTimes(i,1) & ...
                    all_times_all <= fakePreIcTimes(i,2));
                fakePreIcSpikes = [fakePreIcSpikes;spike_idx];
            end
            
            fakeInterIcSpikes = [];
            for i = 1:size(fakeInterIcTimes,1)
                spike_idx = find(all_times_all >= fakeInterIcTimes(i,1) & ...
                    all_times_all <= fakeInterIcTimes(i,2));
                fakeInterIcSpikes = [fakeInterIcSpikes;spike_idx];
            end
            
            % Now get the cluster distributions
            fakePreIcClust = idx(fakePreIcSpikes);
            fakeInterIcClust = idx(fakeInterIcSpikes);
            
            % Get the chi squared
            [~,chi2_boot(ib)] = crosstab([ones(length(fakePreIcClust),1);...
                2*ones(length(fakeInterIcClust),1)],...
                [fakePreIcClust;fakeInterIcClust]);  
        end
        
        % Do it once for real
        realPreIcSpikes = [];
        for i = 1:size(preIcTimesQI,1)
            spike_idx = find(all_times_all >= preIcTimesQI(i,1) & ...
                all_times_all <= preIcTimesQI(i,2));
            realPreIcSpikes = [realPreIcSpikes;spike_idx];
        end

        realInterIcSpikes = [];
        for i = 1:size(interIcTimesQI,1)
            spike_idx = find(all_times_all >= interIcTimesQI(i,1) & ...
                all_times_all <= interIcTimesQI(i,2));
            realInterIcSpikes = [realInterIcSpikes;spike_idx];
        end
        
        realPreIcClust = idx(realPreIcSpikes);
        realInterIcClust = idx(realInterIcSpikes);
        
         % Get the chi squared
        [tbl_2,chi2_real] = crosstab([ones(length(realPreIcClust),1);...
            2*ones(length(realInterIcClust),1)],...
            [realPreIcClust;realInterIcClust]); 
        
        sorted_boot = sort(chi2_boot);
        diff_s_chi = sorted_boot-chi2_real;
        allLarger = find(diff_s_chi>0);
        
        % The p value is the percentage larger 
        if isempty(allLarger) == 1
            p_2 = 1/(nboot+1);
        else
            firstLarger = allLarger(1);
            p_2 = (nboot-firstLarger+1)/(nboot+1);
        end
        
        if k == length(bad_cluster) + 1
            p_2 = 1;
            fprintf('Only one cluster, defining p-value to be 1\n');
        end
        
        chi_tables_plot{whichPt} = tbl_2;
        p_plot(whichPt) = p_2;
        
        % Plot the result of the permuation test
        if doPermPlot == 1
            figure 
            
            subplot(1,2,1)
            hist(sorted_boot)
            hold on
            plot([chi2_real chi2_real],ylim);
            text(mean(xlim),mean(ylim),sprintf('%1.1e',p_2));
            
            % Plot a chi2 distribution
            subplot(1,2,2)
            df = (size(tbl_2,1)-1)*(size(tbl_2,2)-1); % degrees of freedom
            fake_chi_2 = chi2rnd(df,nboot,1);
            hist(fake_chi_2);

            pause
            close(gcf)
        end
        
        
        %% Do post-ictal
        % Concatenate all continuous blocks of time. Very commonly a
        % preictal block will be immediately be preceded by an interictal
        % block. I want to smush those together
        allTimes = [interIcTimesQI;preIcTimesQI;postIcTimesQI];
        allTimes = sortrows(allTimes,1);
        allowableBlocks = allTimes(1,:);
        for i = 2:size(allTimes,1)
            % If the first time here is very close to the last time of the
            % last allowable block, concatenate them
            if abs(allTimes(i,1) - allowableBlocks(end,2)) < 60 % 60 seconds apart
                allowableBlocks(end,2) = allTimes(i,2);
                
            % If not, add this row as a next chunk of allowable time
            else
                allowableBlocks = [allowableBlocks;allTimes(i,:)];
            end
        end
        
        % I now have an array of all allowable chunks of time. My plan will
        % be for each permutation to pick an integer going from 1 to
        % allTotalTime. I then find which cumsum it is in and that tells me
        % which block to try to stick it in.
        totalTimeInChunks = diff(allowableBlocks,1,2);
        allTotalTime = round(sum(totalTimeInChunks));
        cumSumTimes = cumsum(totalTimeInChunks);
        
        
        % Sort postIcTimes by duration
        postIcTimesDuration = diff(postIcTimesQI,1,2);
        [postIcTimesDurationSorted,I] = sort(postIcTimesDuration,'descend');
        postIcTimesSorted = postIcTimesQI(I,:);
        
        nboot = 1e3;
        chi2_boot = zeros(nboot,1);
        
        for ib = 1:nboot
            
            if mod(ib,100) == 0
                fprintf('Doing %d of %d\n',ib,nboot);
            end
            
            fakePostIcTimes = zeros(size(postIcTimesSorted));
            
            % which Post-ictal time am I faking?
            for j = 1:size(fakePostIcTimes,1)
                iter = 0;
                while 1
                    iter = iter + 1;
                    if iter == 1e3
                        error('WHat\n');
                    end
                        
                    t = randi(allTotalTime);

                    % Find which cumsum t fits into. It fits into the last
                    % non-negative one
                    difft = find(t-cumSumTimes>0);
                    if isempty(difft) == 1
                        whichBlock = 1;
                    else
                        whichBlock = min(difft(end)+1,size(allowableBlocks,1));
                    end
                    allowableTimes = allowableBlocks(whichBlock,:);
                    blockDur = diff(allowableTimes);

                    % Now I know which allowable block to stick it in. Pick
                    % a random time in that block, far enough from the end.
                    currDur = postIcTimesDurationSorted(j); % Duration of current Postictal block
                    timeLeft = blockDur -  currDur;
                    if timeLeft < 0
                        continue
                    end
                    t2 = randi(round(timeLeft+1));
                    
                    % The candidate fake Post-ictal time
                    fakePostIcTime = [allowableTimes(1) + t2,...
                        allowableTimes(1) + t2 + currDur];
                    
                    % Throw it out and try again if for some reason it goes
                    % beyond the end of the block (it shouldn't)
                    if fakePostIcTime(2) > allowableTimes(2)+2
                        continue
                    end
                    
                    intersect = 0;
                    
                    % Throw it out and try again if it intersects any of
                    % the previous fake Post-ictal times
                    for kk = 1:j-1
                        oldRange = fakePostIcTimes(kk,:);
                        
                        % If intersect
                        if doTimeRangesIntersect(fakePostIcTime,oldRange) == 1
                            intersect = 1;
                            break
                        end
                        
                    end
                    
                    if intersect == 1, continue; end
                    
                    % If it makes to here, it is acceptable
                    fakePostIcTimes(j,:) = fakePostIcTime;
                    break;
                
                end
                
            end
            
            % Resort the fake PostIctal times
            fakePostIcTimes = sortrows(fakePostIcTimes);
            
            % Do a check for intersections
            intersection = 0;
            
            for i = 1:size(fakePostIcTimes)
                for j = i+1:size(fakePostIcTimes)
                    if doTimeRangesIntersect(fakePostIcTimes(i,:),...
                            fakePostIcTimes(j,:)) == 1
                        intersection = 1;
                    end
                end
            end
            
            if intersection == 1
                error('What\n');
            end
           
            % Now I need to build the fake other times
            fakeOtherIcTimes = [];
            
            % Start with the first block and the first post-ictal time
            startTime = allowableBlocks(1,1);
            whichBlock = 1;
            whichPostIc = 1;
            
            while 1
                currRange = [startTime allowableBlocks(whichBlock,2)];
                if doTimeRangesIntersect(currRange,fakePostIcTimes(whichPostIc,:)) == 1
                    fakeOtherIcTimes = [fakeOtherIcTimes;...
                        startTime fakePostIcTimes(whichPostIc,1)];
                    startTime = fakePostIcTimes(whichPostIc,2);
                    whichPostIc = whichPostIc + 1;
                else
                    whichPostIc = whichPostIc + 1;
                end
                
                % If I've reached the end of the post-ictal blocks
                if whichPostIc == size(fakePostIcTimes,1) + 1
                    
                    % Add the remaining times
                    fakeOtherIcTimes = [fakeOtherIcTimes;...
                            startTime allowableBlocks(whichBlock,2)];
                    
                    % Increase the block and start again
                    whichBlock = whichBlock + 1;
                    if whichBlock == size(allowableBlocks,1) + 1
                        break
                    end
                    startTime = allowableBlocks(whichBlock,1);
                    whichPostIc = 1;
                        
                end 
            end
            
            % Do a check for intersections
            intersection = 0;
            
            for i = 1:size(fakeOtherIcTimes)
                for j = i+1:size(fakeOtherIcTimes)
                    if doTimeRangesIntersect(fakeOtherIcTimes(i,:),...
                            fakeOtherIcTimes(j,:)) == 1
                        intersection = 1;
                    end
                end
            end
            
            if intersection == 1
                error('What\n');
            end
            
            % Do a final check for any intersections
            intersection = 0;
            
            for i = 1:size(fakeOtherIcTimes)
                for j = 1:size(fakePostIcTimes)
                    if doTimeRangesIntersectForgiving(fakeOtherIcTimes(i,:),...
                            fakePostIcTimes(j,:)) == 1
                        intersection = 1;
                    end
                end
            end
            
            if intersection == 1
                error('What\n');
            end
            
            % Do a final check to make sure times add up
            totalTimeFake = sum(diff(fakeOtherIcTimes,1,2)) + ...
                sum(diff(fakePostIcTimes,1,2));
            if abs(totalTimeFake-allTotalTime) > 2
                error('What\n');
            end
            
            
            % Plot the new times
            if 1 == 0
                
                f=figure;
                
                
                for j = 1:size(fakePostIcTimes,1)
                   area([fakePostIcTimes(j,1) fakePostIcTimes(j,2)]/3600,[1 1],'FaceColor','g');
                   hold on
                end
                
                
                for j = 1:size(fakeOtherIcTimes,1)
                   area([fakeOtherIcTimes(j,1) fakeOtherIcTimes(j,2)]/3600,[1 1],'FaceColor','r');
                   hold on
                end
                
                yl = ylim;
                for j = 1:size(szTimes,1)
                   plot([szTimes(j,1) szTimes(j,1)]/3600,yl,'k--','LineWidth',5);
                end
                pause
                close(f);
                
            end
            
            % Now I have the correct times, let me get the spikes in the
            % times
            fakePostIcSpikes = [];
            for i = 1:size(fakePostIcTimes,1)
                spike_idx = find(all_times_all >= fakePostIcTimes(i,1) & ...
                    all_times_all <= fakePostIcTimes(i,2));
                fakePostIcSpikes = [fakePostIcSpikes;spike_idx];
            end
            
            fakeOtherIcSpikes = [];
            for i = 1:size(fakeOtherIcTimes,1)
                spike_idx = find(all_times_all >= fakeOtherIcTimes(i,1) & ...
                    all_times_all <= fakeOtherIcTimes(i,2));
                fakeOtherIcSpikes = [fakeOtherIcSpikes;spike_idx];
            end
            
            % Now get the cluster distributions
            fakePostIcClust = idx(fakePostIcSpikes);
            fakeOtherIcClust = idx(fakeOtherIcSpikes);
            
            % Get the chi squared
            [~,chi2_boot(ib)] = crosstab([ones(length(fakePostIcClust),1);...
                2*ones(length(fakeOtherIcClust),1)],...
                [fakePostIcClust;fakeOtherIcClust]);  
        end
        
        % Do it once for real
        realPostIcSpikes = [];
        for i = 1:size(postIcTimesQI,1)
            spike_idx = find(all_times_all >= postIcTimesQI(i,1) & ...
                all_times_all <= postIcTimesQI(i,2));
            realPostIcSpikes = [realPostIcSpikes;spike_idx];
        end

        realOtherIcSpikes = [];
        otherTimesQI = [preIcTimesQI;interIcTimesQI];
        for i = 1:size(otherTimesQI,1)
            spike_idx = find(all_times_all >= otherTimesQI(i,1) & ...
                all_times_all <= otherTimesQI(i,2));
            realOtherIcSpikes = [realOtherIcSpikes;spike_idx];
        end
        
        realPostIcClust = idx(realPostIcSpikes);
        realOtherIcClust = idx(realOtherIcSpikes);
        
         % Get the chi squared
        [tbl_2,chi2_real] = crosstab([ones(length(realPostIcClust),1);...
            2*ones(length(realOtherIcClust),1)],...
            [realPostIcClust;realOtherIcClust]); 
        
        sorted_boot = sort(chi2_boot);
        diff_s_chi = sorted_boot-chi2_real;
        allLarger = find(diff_s_chi>0);
        
        % The p value is the percentage larger 
        if isempty(allLarger) == 1
            p_2 = 1/(nboot+1);
        else
            firstLarger = allLarger(1);
            p_2 = (nboot-firstLarger+1)/(nboot+1);
        end
        
        if k == length(bad_cluster) + 1
            p_2 = 1;
            fprintf('Only one cluster, defining p-value to be 1\n');
        end
        
        chi_tables_postIc{whichPt} = tbl_2;
        p_postIc(whichPt) = p_2;
        
        % Plot the result of the permuation test
        if doPermPlot == 1
            figure 
            
            subplot(1,2,1)
            hist(sorted_boot)
            hold on
            plot([chi2_real chi2_real],ylim);
            text(mean(xlim),mean(ylim),sprintf('%1.1e',p_2));
            
            % Plot a chi2 distribution
            subplot(1,2,2)
            df = (size(tbl_2,1)-1)*(size(tbl_2,2)-1); % degrees of freedom
            fake_chi_2 = chi2rnd(df,nboot,1);
            hist(fake_chi_2);

            pause
            close(gcf)
        end
        
    end
end

if intericTime == 4
    save([destFolder,'stats4TimePerm.mat'],'stats');
elseif intericTime == 1
    save([destFolder,'stats1TimePerm.mat'],'stats');
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
        print(gcf,[destFolder,'clustBarTimePerm',sprintf('%d_hours',intericTime)],'-depsc');
        eps2pdf([destFolder,'clustBarTimePerm',sprintf('%d_hours',intericTime),'.eps'])
        end

        %% Post-ic
        
        all_p = [];
        for whichPt = whichPts
           all_p = [all_p;p_postIc(whichPt)]; 
        end

        X_2 = -2 * sum(log(all_p));
        sum_p = 1-chi2cdf(X_2,2*length(all_p));

        fprintf('The group p value for post-ictal vs interictal cluster is %1.1e\n',sum_p);
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
        print(gcf,[destFolder,'clustPostTimePerm',sprintf('%d_hours',intericTime)],'-depsc');
        eps2pdf([destFolder,'clustPostTimePerm',sprintf('%d_hours',intericTime),'.eps'])


    end
   

end


end