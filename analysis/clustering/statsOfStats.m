function statsOfStats(pt,stats)

p_pre_all = [];
p_post_all = [];
p_dist_all = [];
p_hour_all = [];
whichPts = [];
diff_dist_all = [];
outcome_all = [];
temp_lobe_all = [];

for whichPt = 1:length(stats)
    if isempty(stats(whichPt).hour), continue; end
    
    whichPts = [whichPts,whichPt];
    
    
    % Get stats on if it changes from hour to hour
    p_hour = stats(whichPt).hour.p;
    p_hour_all = [p_hour_all,p_hour];
    
    
    % Get stats on if it changes in pre-ictal period
    p_pre = stats(whichPt).preIc.p;
    % P can't be lower than the number of permutations
    if p_pre == 0
        p_pre = 1/1000;
    end
    p_pre_all = [p_pre_all,p_pre];
    
    % Get stats on if it changes in post ictal period
    p_post = stats(whichPt).postIc.p;
    if p_post == 0
        p_post = 1/1000;
    end
    p_post_all = [p_post_all,p_post];
    
    % Get stats on distance change
    p_dist = stats(whichPt).dist.p;
    if p_dist == 0
        p_dist = 1/1000;
    end
    if isempty(p_dist)==1, p_dist = nan; end
    
    if 1==1%p_post < 0.05/20
        p_dist_all = [p_dist_all,p_dist];
        dist_diff = stats(whichPt).dist.dist_diff;
        if isempty(dist_diff) == 1, dist_diff = nan; end
        diff_dist_all = [diff_dist_all,dist_diff];
    end
    
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

%% Pre-ic change
X_2 = -2 * sum(log(p_post_all));
sum_p = 1-chi2cdf(X_2,2*length(p_post_all));
fprintf('%d out of %d patients showed a significant pre-ictal change in cluster distribution.\n',...
    sum(p_post_all < 0.05/length(whichPts)),length(whichPts))
fprintf('The group p value for change over time is %1.1e\n',sum_p);


%% Distance change


%% Clinical correlations

% correlate whether there is an hour-to-hour change with clinical outcome
hourChange = p_hour_all < 0.05/length(p_hour_all);
[p_hour_outcome,info_hour_outcome] = correlateClinically(hourChange,outcome_all,'bin','num',1);

% correlate whether there is an hour-to-hour change with temporal/non
% temporal lobe
[p_hour_lobe,info_hour_lobe] = correlateClinically(hourChange,temp_lobe_all,'bin','bin',1);

end