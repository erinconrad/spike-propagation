function statsOfStats(pt,stats)

%% Initialize stats variables for all patients
names = {};
p_pre_all = [];
p_post_all = [];
p_dist_all = [];
p_hour_all = [];
whichPts = [];
soz_dist_pre_all = [];
soz_dist_inter_all = [];
soz_dist_post_all = [];
soz_dist_other_all = [];
SD_pre_all = [];
SD_post_all = [];
SD_other_all = [];
SD_inter_all = [];
outcome_all = [];
temp_lobe_all = [];

%% Loop through patients and grab stats info
for whichPt = 1:length(stats)
    if isempty(stats(whichPt).hour), continue; end
    
    whichPts = [whichPts,whichPt];
    names = [names,pt(whichPt).name];
    
    % Get stats on if it changes from hour to hour
    p_hour = stats(whichPt).hour.p;
    p_hour_all = [p_hour_all,p_hour];
    
    % Get stats on if it changes in pre-ictal period
    p_pre = stats(whichPt).preIc.p;
    % P can't be lower than the number of permutations
    if p_pre == 0
        p_pre = 1/1001;
    end
    p_pre_all = [p_pre_all,p_pre];
    
    % Get stats on if it changes in post ictal period
    p_post = stats(whichPt).postIc.p;
    if p_post == 0
        p_post = 1/1001;
    end
    p_post_all = [p_post_all,p_post];
    
    % Distance from SOZ
    pre_dist = stats(whichPt).soz.pre.pre_dist;
    inter_dist = stats(whichPt).soz.pre.inter_dist;
    post_dist = stats(whichPt).soz.post.post_dist;
    other_dist = stats(whichPt).soz.post.other_dist;
    
    soz_dist_pre_all = [soz_dist_pre_all,pre_dist];
    soz_dist_inter_all = [soz_dist_inter_all,inter_dist];
    soz_dist_post_all = [soz_dist_post_all,post_dist];
    soz_dist_other_all = [soz_dist_other_all,other_dist];
    
    
    % Standard distance (spatial dispersion)
    SD_pre = stats(whichPt).dispersion.pre.SD_pre;
    SD_inter = stats(whichPt).dispersion.pre.SD_inter;
    SD_post = stats(whichPt).dispersion.post.SD_post;
    SD_other = stats(whichPt).dispersion.post.SD_other;
    
    SD_pre_all = [SD_pre_all,SD_pre];
    SD_inter_all = [SD_inter_all,SD_inter];
    SD_post_all = [SD_post_all,SD_post];
    SD_other_all = [SD_other_all,SD_other];
    
    
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
 
end

%% Hour change
X_2 = -2 * sum(log(p_hour_all));
sum_p = 1-chi2cdf(X_2,2*length(p_hour_all));
fprintf('%d out of %d patients showed a significant change in cluster distribution over time.\n',...
    sum(p_hour_all < 0.05/length(whichPts)),length(whichPts))
fprintf('The group p value for change over time is %1.1e\n',sum_p);


%% Pre-ic change
X_2 = -2 * sum(log(p_pre_all));
sum_p = 1-chi2cdf(X_2,2*length(p_pre_all));
fprintf('%d out of %d patients showed a significant pre-ictal change in cluster distribution.\n',...
    sum(p_pre_all < 0.05/length(whichPts)),length(whichPts))
fprintf('The group p value for change over time is %1.3f\n',sum_p);

%% Post-ic change
X_2 = -2 * sum(log(p_post_all));
sum_p = 1-chi2cdf(X_2,2*length(p_post_all));
fprintf('%d out of %d patients showed a significant post-ictal change in cluster distribution.\n',...
    sum(p_post_all < 0.05/length(whichPts)),length(whichPts))
fprintf('The group p value for change over time is %1.1e\n',sum_p);


%% Change in distance from SOZ
[~,p] = ttest2(soz_dist_pre_all,soz_dist_inter_all);
fprintf(['The p-value for change in distance from SOZ in the'...
    'pre-ictal period is:\n%1.2e\n'],p);

[~,p] = ttest2(soz_dist_post_all,soz_dist_other_all);
fprintf(['The p-value for change in distance from SOZ in the'...
    'post-ictal period is:\n%1.2e\n'],p);

%% Change in SD
[~,p] = ttest2(SD_pre_all,SD_inter_all);
fprintf(['The p-value for the change in standard distance in the'...
    'pre-ictal period is:\n%1.2e\n'],p);

[~,p] = ttest2(SD_post_all,SD_other_all);
fprintf(['The p-value for the change in standard distance in the'...
    'post-ictal period is:\n%1.2e\n'],p);

%% Clinical correlations

% correlate whether there is an hour-to-hour change with clinical outcome
hourChange = p_hour_all < 0.05/length(p_hour_all);
[p_hour_outcome,info_hour_outcome] = correlateClinically(hourChange,outcome_all,'bin','num',0);
fprintf(['The p-value for Wilcoxon rank sum comparing outcome between\n'...
    'patients with hour-to-hour change and those without is:/n'...
    'p = %1.2e\n'],p_hour_outcome);

% correlate whether there is an hour-to-hour change with temporal/non
% temporal lobe
[p_hour_lobe,info_hour_lobe] = correlateClinically(hourChange,temp_lobe_all,'bin','bin',0);
fprintf(['The p-value for chi squared comparing temporal vs non-temporal lobe between\n'...
    'patients with hour-to-hour change and those without is:/n'...
    'p = %1.2e\n'],p_hour_lobe);

end