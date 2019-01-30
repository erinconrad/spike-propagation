function AD_AR(pt,cluster,power,whichPts)

plotInfo = 0;

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
p1 = genpath(scriptFolder);
addpath(p1);
destFolder = [resultsFolder,'alphaDelta/plots/'];
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


allP = [];
allT = [];
outcome_all = [];
temp_lobe_all = [];

all_b_soz = [];
all_p_soz = [];
all_b_SL = [];
all_p_SL = [];
names = {};


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
    
    % outcomes
    outcome = getOutcome(pt,whichPt);
    outcome_all = [outcome_all,outcome];
    
    % SOZ
    szOnsetText = pt(whichPt).clinical.seizureOnset;
    if contains(szOnsetText,'TL') == 1
        tempLobe = 1;
    else
        tempLobe = 0;
    end
    temp_lobe_all = [temp_lobe_all,tempLobe];

    
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
    
    % Get most popular cluster
    popular = mode(idx);
    
   
    
    %% Get the sequence lengths of all spikes
    [seq_lengths,seq_times] = getSeqDist(pt,cluster,whichPt);
    
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
    
    % Run through bin times and get proportion of spikes in most popular
    % cluster for that bin.
    for i = 1:size(bin_times,1)
        
        % Cluster identities of spikes in between those times
        whichClust = idx(all_times_all > bin_times(i,1) & ...
            all_times_all < bin_times(i,2));
        
        prop_pop(i) = sum(whichClust == popular)/length(whichClust);
        
        
        % Get spike locs
        whichSpikes = find(all_times_all > bin_times(i,1) & ...
            all_times_all < bin_times(i,2));
        locs_bin(i,:) = mean(all_locs(whichSpikes,:),1);
        num_spikes(i) = length(whichSpikes);
        
        
        
        % Get mean sequence length in these times
        SL_bin(i) = mean(seq_lengths(seq_times > bin_times(i,1) & ...
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
    
    %% Remove nans for distance from SOZ analysis
    nan_times = find(isnan(soz_dist_bin));
    times = old_times;
    mean_ad = old_mean_ad;
    
    times(nan_times) = [];
    mean_ad(nan_times) = [];
    soz_dist_bin(nan_times) = [];
    
    %% Do initial regression
    Y = soz_dist_bin;
    
    % X is the predictor. The first component of X is the alpha-delta
    % ratio, the predictor I am interested in. The next component is a
    % constant error term. The third component is just what time it is,
    % reflecting a linear trend with time. The last is the categorical
    % variable representing the hour of the day, reflecting a cyclical Q24
    % hour trend.
    %X = [mean_ad ones(size(mean_ad)) times cat_hours];
    X = [mean_ad ones(size(mean_ad))];
    
    [b,~,resid,~,stats] = regress(Y, X);
    
    %% Get initial guess for autocorrelation term
    r = corr(resid(1:end-1),resid(2:end));  
    
    %% Define anonymous function representing new fit including autocorrelation
    f = @(c,x) [Y(1); c(1)*Y(1:end-1) + (x(2:end,:)- c(1)*x(1:end-1,:))*c(2:end)];
    [c,~,~,CovB,~,~] = nlinfit(X,Y,f,[r;b]);
    mdl = fitnlm(X,Y,f,[r;b]);
    p = mdl.Coefficients.pValue(2);
    r2 = mdl.Rsquared.Adjusted;
    
    % c(1) is the autocorrelation term
    % c(2) is b(1) which is the linear correlation term
    % c(3) is b(2) which is the constant error term
    
    %p = linhyptest(c(3),CovB(3,3));
    
    [~,p_rank] = corr(mean_ad,soz_dist_bin,'Type','Spearman');
    
    if plotInfo == 1
        
        fprintf(['For %s, the rank p-value is %1.1e,\n'...
            'and the non-linear is %1.1e.\n\n\n\n\n'],...
            pt(whichPt).name,p_rank,p);
        
        %Info about new residuals
        figure
        subplot(1,3,1)
        u = Y - f(c,X);
        plot(u);
        title('Residuals');
        
        subplot(1,3,2)
        plot(times,soz_dist_bin,'b');
        hold on
        plot(times,X*b,'r');
        legend('Real dist from SOZ','original model');
        
        subplot(1,3,3)
        fakeY = f(c,X);
        plot(times,soz_dist_bin,'b');
        hold on
        plot(times,fakeY,'g');
        legend('Real dist from SOZ','non linear model');
        pause
        close(gcf)
    end
    
    all_b_soz = [all_b_soz;c(2)];
    all_p_soz = [all_p_soz;p];
    
    
    
    %% Now do SL
    
    nan_times = find(isnan(SL_bin));
    times = old_times;
    mean_ad = old_mean_ad;
    
    times(nan_times) = [];
    mean_ad(nan_times) = [];
    SL_bin(nan_times) = [];
    
    %% Do initial regression
    Y = SL_bin;
    
    % X is the predictor. The first component of X is the alpha-delta
    % ratio, the predictor I am interested in. The next component is a
    % constant error term. The third component is just what time it is,
    % reflecting a linear trend with time. The last is the categorical
    % variable representing the hour of the day, reflecting a cyclical Q24
    % hour trend.
    %X = [mean_ad ones(size(mean_ad)) times cat_hours];
    X = [mean_ad ones(size(mean_ad))];
    
    [b,~,resid,~,stats] = regress(Y, X);
    
    %% Get initial guess for autocorrelation term
    r = corr(resid(1:end-1),resid(2:end));  
    
    %% Define anonymous function representing new fit including autocorrelation
    f = @(c,x) [Y(1); c(1)*Y(1:end-1) + (x(2:end,:)- c(1)*x(1:end-1,:))*c(2:end)];
    [c,~,~,CovB,~,~] = nlinfit(X,Y,f,[r;b]);
    mdl = fitnlm(X,Y,f,[r;b]);
    p = mdl.Coefficients.pValue(2);
    r2 = mdl.Rsquared.Adjusted;
    t = mdl.Coefficients.tStat(2);
    
    % c(1) is the autocorrelation term
    % c(2) is b(1) which is the linear correlation term
    % c(3) is b(2) which is the constant error term
    
    %p = linhyptest(c(3),CovB(3,3));
    
    [~,p_rank] = corr(mean_ad,SL_bin,'Type','Spearman');
    
    if plotInfo == 1
        
        fprintf(['For %s, the rank p-value is %1.1e,\n'...
            'and the non-linear is %1.1e.\n\n\n\n\n'],...
            pt(whichPt).name,p_rank,p);
        
        %Info about new residuals
        figure
        subplot(1,3,1)
        u = Y - f(c,X);
        plot(u);
        title('Residuals');
        
        subplot(1,3,2)
        plot(times,SL_bin,'b');
        hold on
        plot(times,X*b,'r');
        legend('Real sequence length','original model');
        
        subplot(1,3,3)
        fakeY = f(c,X);
        plot(times,SL_bin,'b');
        hold on
        plot(times,fakeY,'g');
        legend('Real sequence length','non linear model');
        pause
        close(gcf)
    end
    
    all_b_SL = [all_b_SL;c(2)];
    all_p_SL = [all_p_SL;p];
   
    
    
    %% DO PROP-POP ANALYSIS
    
    if length(clusters) == 1
        fprintf('One cluster for %s, skipping\n',pt(whichPt).name);
        p = 1;
        allP =[allP;p];
        allT = [allT;nan];
        continue
    end
    
    %% Remove nan time points for prop-pop analysis
    nan_times = find(isnan(prop_pop));
    times = old_times;
    mean_ad = old_mean_ad;
    
    
    times(nan_times) = [];
    mean_ad(nan_times) = [];
    prop_pop(nan_times) = [];
    
    
    %% Do initial regression incorporating linear trend and Q24 hour trend
    
    % Get hours
    hours = floor(mod(times,24*3600)/3600) + 1;
    
    % Get hours as a categorical variable
    cat_hours = dummyvar(hours);
    
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
    
    % Do the regression
    [b,bint,resid,~,stats] = regress(Y, X);
    
    
    %Y_fake = X*b;
    %p = dwtest(resid,X);
    
    %{
    figure
    plot(resid)
    pause
    close(gcf)
    %}
    
    %% Now correct for autocorrelation
    
    % Get initial guess for autocorrelation
    r = corr(resid(1:end-1),resid(2:end));  
    
    % Define anonymous function representing new fit including autocorrelation
    f = @(c,x) [Y(1); c(1)*Y(1:end-1) + (x(2:end,:)- c(1)*x(1:end-1,:))*c(2:end)];
    
    % Do the new model 
    [c,~,~,CovB,~,~] = nlinfit(X,Y,f,[r;b]);
    mdl = fitnlm(X,Y,f,[r;b]);
    
    % c(1) is the autocorrelation term
    % c(2) is b(1) which is the linear correlation term
    % c(3) is b(2) which is the constant error term
    % c(4) is b(3) which is linear trend
    % c(5-end) is b(4-end) which are the seasonal components
    
   % p = linhyptest(c(3),CovB(3,3));
    
    p = mdl.Coefficients.pValue(2);
    t = mdl.Coefficients.tStat(2);
    r2 = mdl.Rsquared.Adjusted;
    
    allP = [allP;p];
    allT = [allT;t];
    
    [~,p_rank] = corr(mean_ad,prop_pop,'Type','Spearman');
    
    if plotInfo == 1
        
        
        fprintf(['For %s, the rank p-value is %1.1e,\n'...
            'and the non-linear is %1.1e.\n']...
            ,pt(whichPt).name,p_rank,p);
        
        %Info about new residuals
        figure
        subplot(1,3,1)
        u = Y - f(c,X);
        plot(u);
        title('Residuals');
        
        subplot(1,3,2)
        plot(times,prop_pop,'b');
        hold on
        plot(times,X*b,'r');
        legend('Real prop-pop','original model');
        
        subplot(1,3,3)
        plot(times,prop_pop,'b');
        hold on
        fakeY = f(c,X);
        plot(times,fakeY)
        legend('Real prop-pop','non-linear model');
        
        pause
        close(gcf)
    end
    
    
    
    
    
    
    
    
    
    
end


X_2 = -2 * sum(log(allP));
sum_p = 1-chi2cdf(X_2,2*length(allP));
fprintf(['There are %d with significant correlation.\nThe combined p-value'...
    'is %1.1e.\n'],sum(allP<0.05/length(allP)),sum_p);

%% Make table of p-values and t-statistics
p_text = getPText(allP);
allT_text = num2str(allT,3);
table(char(names),allT_text,char(p_text))

%% Test that b significantly different from zero for SOZ
[~,p,ci,stats] = ttest(all_b_soz);
fprintf('P value for SOZ is %1.1e.\n',p);
changePos = find(allP<0.05/length(allP));
all_b_soz_changePos = all_b_soz(changePos);
all_p_soz_changePos = all_p_soz(changePos);
table(names(changePos),all_b_soz_changePos,all_p_soz_changePos)

%% Correlation between change in location and outcome
%
changeLoc = allP < 0.05/length(allP);
[p,info] = correlateClinically(changeLoc,outcome_all,'bin','num',0);
fprintf(['The p-value for Wilcoxon rank sum comparing outcome between\n'...
    'patients with change and those without is:\n'...
    'p = %1.2e\n'],p);


%% Correlation between change in location and temp lobe
[p,info] = correlateClinically(changeLoc,temp_lobe_all,'bin','bin',0);
fprintf(['The p-value for chi squared comparing temporal vs non-temporal lobe between\n'...
    'patients with change and those without is:\n'...
    'p = %1.2e\n'],p);


%% Table of stuff
p_sleep_t = getPText(allP);
T = table(p_sleep_t);


%% Test that b significantly different from zero for SL
[~,p,ci,stats] = ttest(all_b_SL);
fprintf('P value for SL is %1.1e.\n',p);
changePos = find(allP<0.05/length(allP));
all_b_SL_changePos = all_b_SL(changePos);
all_p_SL_changePos = all_p_SL(changePos);
table(names(changePos),all_b_SL_changePos,all_p_SL_changePos)

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

