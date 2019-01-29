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
SL_pre_all = [];
SL_post_all = [];
SL_inter_all = [];
SL_other_all = [];
SL_p_pre_all = [];
SL_p_post_all = [];

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
    if isfield(stats(whichPt),'soz') == 1
        pre_dist = stats(whichPt).soz.pre.pre_dist;
        inter_dist = stats(whichPt).soz.pre.inter_dist;
        post_dist = stats(whichPt).soz.post.post_dist;
        other_dist = stats(whichPt).soz.post.other_dist;

        soz_dist_pre_all = [soz_dist_pre_all,pre_dist];
        soz_dist_inter_all = [soz_dist_inter_all,inter_dist];
        soz_dist_post_all = [soz_dist_post_all,post_dist];
        soz_dist_other_all = [soz_dist_other_all,other_dist];
    end
    
    
    % Standard distance (spatial dispersion)
    if isfield(stats(whichPt),'dispersion') == 1
        SD_pre = stats(whichPt).dispersion.pre.SD_pre;
        SD_inter = stats(whichPt).dispersion.pre.SD_inter;
        SD_post = stats(whichPt).dispersion.post.SD_post;
        SD_other = stats(whichPt).dispersion.post.SD_other;

        SD_pre_all = [SD_pre_all,SD_pre];
        SD_inter_all = [SD_inter_all,SD_inter];
        SD_post_all = [SD_post_all,SD_post];
        SD_other_all = [SD_other_all,SD_other];
    end
    
    
    % Sequence length
    if isfield(stats(whichPt),'SL') == 1
        SL_pre = stats(whichPt).SL.pre.SL_pre;
        SL_inter = stats(whichPt).SL.pre.SL_inter;
        SL_p_pre = stats(whichPt).SL.pre.p;
        
        SL_post = stats(whichPt).SL.post.SL_post;
        SL_other = stats(whichPt).SL.post.SL_other;
        SL_p_post = stats(whichPt).SL.post.p;
        
        SL_pre_all = [SL_pre_all;SL_pre];
        SL_inter_all = [SL_inter_all;SL_inter];
        SL_post_all = [SL_post_all;SL_post];
        SL_other_all = [SL_other_all;SL_other];
        SL_p_pre_all = [SL_p_pre_all;SL_p_pre];
        SL_p_post_all = [SL_p_post_all;SL_p_post];
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

%% Post-ic change
X_2 = -2 * sum(log(p_post_all));
sum_p = 1-chi2cdf(X_2,2*length(p_post_all));
fprintf('%d out of %d patients showed a significant post-ictal change in cluster distribution.\n',...
    sum(p_post_all < 0.05/length(whichPts)),length(whichPts))
fprintf('The group p value for change over time is %1.1e\n',sum_p);


%% Change in distance from SOZ
if isempty(soz_dist_pre_all) == 0
    [~,p] = ttest2(soz_dist_pre_all,soz_dist_inter_all);
    fprintf(['The p-value for change in distance from SOZ in the'...
        'pre-ictal period is:\n%1.2e\n'],p);

    [~,p] = ttest2(soz_dist_post_all,soz_dist_other_all);
    fprintf(['The p-value for change in distance from SOZ in the'...
        'post-ictal period is:\n%1.2e\n'],p);

    [p,tbl_kw,stats_kw] = kruskalwallis([soz_dist_pre_all',soz_dist_inter_all',...
        soz_dist_post_all'],[],'off');
    fprintf(['The p-value for change in distance from SOZ comparing'...
        'pre/post/inter by K-W is:\n%1.2e\n'],p);
end


%% Change in SD
if isempty(SD_pre_all) == 0
    [~,p] = ttest2(SD_pre_all',SD_inter_all');
    fprintf(['The p-value for the change in standard distance in the'...
        'pre-ictal period is:\n%1.2e\n'],p);

    [~,p] = ttest2(SD_post_all',SD_other_all');
    fprintf(['The p-value for the change in standard distance in the'...
        'post-ictal period is:\n%1.2e\n'],p);
end

%% Clinical correlations

% correlate whether there is an hour-to-hour change with clinical outcome
hourChange = p_hour_all < 0.05/length(p_hour_all);
[p_hour_outcome,info_hour_outcome] = correlateClinically(hourChange,outcome_all,'bin','num',0);
fprintf(['The p-value for Wilcoxon rank sum comparing outcome between\n'...
    'patients with hour-to-hour change and those without is:\n'...
    'p = %1.2e\n'],p_hour_outcome);

% correlate whether there is an hour-to-hour change with temporal/non
% temporal lobe
[p_hour_lobe,info_hour_lobe] = correlateClinically(hourChange,temp_lobe_all,'bin','bin',0);
fprintf(['The p-value for chi squared comparing temporal vs non-temporal lobe between\n'...
    'patients with hour-to-hour change and those without is:\n'...
    'p = %1.2e\n'],p_hour_lobe);

% Correlate clinical outcome with post-ictal change
postIcChange = p_post_all < 0.05/length(whichPts);
[p_post_outcome,info_post_outcome] = correlateClinically(postIcChange,outcome_all,'bin','num',0);
fprintf(['The p-value for Wilcoxon rank sum comparing outcome between\n'...
    'patients with post-ictal change and those without is:\n'...
    'p = %1.2e\n'],p_post_outcome);
    
% correlate post-ictal change with temporal lobe onset
[p_post_lobe,info_post_lobe] = correlateClinically(postIcChange,temp_lobe_all,'bin','bin',0);
fprintf(['The p-value for chi squared comparing temporal vs non-temporal lobe between\n'...
    'patients with post-ictal change and those without is:\n'...
    'p = %1.2e\n'],p_post_lobe);


%% Table of stuff
p_hour_all_t = getPText(p_hour_all);
p_pre_all_t = getPText(p_pre_all);
p_post_all_t = getPText(p_post_all);

T = table(names',p_hour_all_t,p_pre_all_t,p_post_all_t);



%% SL

if isempty(SL_pre_all) == 0
    [~,p] = ttest2(SL_pre_all,SL_inter_all);
    fprintf(['The p-value for change in SL in the'...
        'pre-ictal period is:\n%1.2e\n'],p);

    [~,p] = ttest2(SL_post_all,SL_other_all);
    fprintf(['The p-value for change in SL in the'...
        'post-ictal period is:\n%1.2e\n'],p);

    
    fprintf('%d of %d patients had a significant pre-ictal change in SL.\n',...
        sum(SL_p_pre_all < 0.05/length(SL_p_pre_all)),length(SL_p_pre_all));
    

    fprintf('%d of %d patients had a significant post-ictal change in SL.\n',...
        sum(SL_p_post_all < 0.05/length(SL_p_post_all)),length(SL_p_post_all));
    
end

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

