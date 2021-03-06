function statsOfStats(pt,stats,cluster)

%{

Does stats on the stats.mat structure that I get from running CNewStats
(cluster statistics)

%}

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
soz_dist_post_p = [];
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
hour_chi2_all = [];
pre_chi2_all = [];
hour_dof_all = [];
pre_dof_all = [];
post_dof_all = [];
post_chi2_all = [];
MTS_all = [];
no_depths_all = [];
age_all = [];

%% Loop through patients and grab stats info
for whichPt = 1:length(stats)
    if isempty(stats(whichPt).hour), continue; end
    
    whichPts = [whichPts,whichPt];
    names = [names,pt(whichPt).name];
    
    % Get stats on if it changes from hour to hour
    p_hour = stats(whichPt).hour.p;
    if p_hour == 0
        p_hour = 1/1001;
    end
    p_hour_all = [p_hour_all,p_hour];
    hour_chi2_all = [hour_chi2_all,stats(whichPt).hour.chi2];
    hour_dof_all = [hour_dof_all,...
        (size(stats(whichPt).hour.tbl,1)-1)*(size(stats(whichPt).hour.tbl,2)-1)];
    
    % Get stats on if it changes in pre-ictal period
    p_pre = stats(whichPt).preIc.p;
    % P can't be lower than the number of permutations
    if p_pre == 0
        p_pre = 1/1001;
    end
    p_pre_all = [p_pre_all,p_pre];
    pre_chi2_all = [pre_chi2_all,stats(whichPt).preIc.chi2];
    pre_dof_all = [pre_dof_all,...
        (size(stats(whichPt).preIc.tbl,1)-1)*(size(stats(whichPt).preIc.tbl,2)-1)];

    
    % Get stats on if it changes in post ictal period
    p_post = stats(whichPt).postIc.p;
    if p_post == 0
        p_post = 1/1001;
    end
    p_post_all = [p_post_all,p_post];
    post_chi2_all = [post_chi2_all,stats(whichPt).postIc.chi2];
    post_dof_all = [post_dof_all,...
        (size(stats(whichPt).postIc.tbl,1)-1)*(size(stats(whichPt).postIc.tbl,2)-1)];
    
    % Distance from SOZ
    if isfield(stats(whichPt),'soz') == 1
        pre_dist = stats(whichPt).soz.pre.pre_dist;
        inter_dist = stats(whichPt).soz.pre.inter_dist;
        post_dist = stats(whichPt).soz.post.post_dist;
        other_dist = stats(whichPt).soz.post.other_dist;
        post_p = stats(whichPt).soz.post.p;

        soz_dist_pre_all = [soz_dist_pre_all,pre_dist];
        soz_dist_inter_all = [soz_dist_inter_all,inter_dist];
        soz_dist_post_all = [soz_dist_post_all,post_dist];
        soz_dist_other_all = [soz_dist_other_all,other_dist];
        soz_dist_post_p = [soz_dist_post_p,post_p];
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
    
    % Age
    age = pt(whichPt).clinical.ageSurgery;
    if contains(age,'-') == 1
        a = regexp(age,'-');
        num1 = str2num(age(1:a-1));
        num2 = str2num(age(a+1:end));
        age = mean([num1,num2]);
    elseif contains(age,'?')
        age = nan;
    else
        age = str2num(age);
    end
    age_all = [age_all;age];
    
    % SOZ
    szOnsetText = pt(whichPt).clinical.seizureOnset;
    if contains(szOnsetText,'TL') == 1
        tempLobe = 1;
    else
        tempLobe = 0;
    end
    temp_lobe_all = [temp_lobe_all,tempLobe];
    
    % Path
    path = pt(whichPt).clinical.pathology;
    if contains(path,'MTS') == 1
        MTS = 1;
    else
        MTS = 0;
    end
    MTS_all = [MTS_all,MTS];
    
    % Electrode type
    [n_grids_all,n_strips_all,n_depths_all] = getElectrodeInfo(pt,cluster,whichPt,0);
    if n_depths_all == 0
        no_depths = 1;
    else
        no_depths = 0;
    end
    no_depths_all = [no_depths_all,no_depths];
 
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
group_pval = fisher_pvalue_meta_analysis(p_post_all);


%% Change in distance from SOZ
if isempty(soz_dist_pre_all) == 0
    [~,p] = ttest2(soz_dist_pre_all,soz_dist_inter_all);
    fprintf(['The p-value for change in distance from SOZ in the'...
        'pre-ictal period is:\n%1.2e\n'],p);

    [~,p] = ttest2(soz_dist_post_all,soz_dist_other_all);
    fprintf(['The p-value for change in distance from SOZ in the'...
        'post-ictal period is:\n%1.2e\n'],p);

    %{
    [p,tbl_kw,stats_kw] = kruskalwallis([soz_dist_pre_all',soz_dist_inter_all',...
        soz_dist_post_all'],[],'off');
    fprintf(['The p-value for change in distance from SOZ comparing'...
        'pre/post/inter by K-W is:\n%1.2e\n'],p);
        %}
        
    [p,tbl,stats] = friedman([soz_dist_pre_all',soz_dist_inter_all',...
        soz_dist_post_all'],1,'off');
    fprintf(['The Friedman test comparing distance from SOZ is:\n'...
        'p = %1.2e, Q = %1.1f\n'],p,tbl{2,5});
    
    fprintf('%d of %d patients had a significant post-ictal change in dist from SOZ.\n',...
        sum(soz_dist_post_p < 0.05/length(soz_dist_post_p)),length(soz_dist_post_p));
    
    
    sig_post_move = find(p_post_all < 0.05/length(whichPts));
    names_post_move = names(sig_post_move)';
    soz_post_move = soz_dist_post_all(sig_post_move);
    soz_other_move = soz_dist_other_all(sig_post_move);
    soz_p_move = soz_dist_post_p(sig_post_move);
    
    table(names_post_move,soz_post_move',soz_other_move',soz_p_move')
    
    
end



%{
%% Change in SD
if isempty(SD_pre_all) == 0
    [~,p] = ttest2(SD_pre_all',SD_inter_all');
    fprintf(['The p-value for the change in standard distance in the'...
        'pre-ictal period is:\n%1.2e\n'],p);

    [~,p] = ttest2(SD_post_all',SD_other_all');
    fprintf(['The p-value for the change in standard distance in the'...
        'post-ictal period is:\n%1.2e\n'],p);
end
%}

%% Clinical correlations

% Overall clinical info
ilae_overall = getILAE(outcome_all);
fprintf('Median ILAE is %1.1f, range %1.1f-%1.1f\n',median(ilae_overall),...
    min(ilae_overall),max(ilae_overall));

fprintf('%d patients had ILAE <= 3\n',sum(ilae_overall<=3));




% correlate whether there is an hour-to-hour change with clinical outcome
hourChange = p_hour_all < 0.05/length(p_hour_all);
[p_hour_outcome,info_hour_outcome] = correlateClinically(hourChange,ilae_overall,'bin','num',0);
[~,~, u_mat] = ranksum_erin(ilae_overall(hourChange==1)',ilae_overall(hourChange==0)');
test_stat = getStandardStats(ilae_overall(hourChange==1)',ilae_overall(hourChange==0)','rs');
fprintf(['The p-value for Wilcoxon rank sum comparing outcome between\n'...
    'patients with hour-to-hour change and those without is:\n'...
    'p = %1.2e\nMatlab U is %1.1f and my U is %1.1f\n'],p_hour_outcome,...
    u_mat,test_stat);

% correlate hour-to-hour change with age
[p_hour_age,info_hour_age] = correlateClinically(hourChange,age_all<=13,'bin','bin',0);
fprintf(['The p-value for chi squared comparing child vs adult between\n'...
    'patients with hour-to-hour change and those without is:\n'...
    'p = %1.2e\n'],p_hour_age);

% Get average ILAE scores of those with a change and those without
ilae_change = getILAE(outcome_all(hourChange == 1));
ilae_no_change = getILAE(outcome_all(hourChange == 0));

fprintf('Median ILAE for change is %1.1f and for no change is %1.1f\n',...
    median(ilae_change),median(ilae_no_change));

% correlate whether there is an hour-to-hour change with temporal/non
% temporal lobe
[p_hour_lobe,info_hour_lobe] = correlateClinically(hourChange,temp_lobe_all,'bin','bin',0);
fprintf(['The p-value for chi squared comparing temporal vs non-temporal lobe between\n'...
    'patients with hour-to-hour change and those without is:\n'...
    'p = %1.2e\n'],p_hour_lobe);

% correlate hour-to-hour change with MTS vs other pathology
[p_hour_MTS,info_hour_MTS] = correlateClinically(hourChange,MTS_all,'bin','bin',0);
fprintf(['The p-value for chi squared comparing MTS vs non-MTS lobe between\n'...
    'patients with hour-to-hour change and those without is:\n'...
    'p = %1.2e, chi2 = %1.2f\n'],p_hour_MTS,info_hour_MTS.chi2);


% correlate hour-to-hour change with no depths vs yes depths
[p_hour_depths,info_hour_depths] = correlateClinically(hourChange,no_depths_all,'bin','bin',0);
fprintf(['The p-value for chi squared comparing depths vs no depths between\n'...
    'patients with hour-to-hour change and those without is:\n'...
    'p = %1.2e, chi2 = %1.2f\n'],p_hour_depths,info_hour_depths.chi2);




% Correlate clinical outcome with post-ictal change
postIcChange = p_post_all < 0.05/length(whichPts);
[p_post_outcome,info_post_outcome] = correlateClinically(postIcChange,ilae_overall,'bin','num',0);
fprintf(['The p-value for Wilcoxon rank sum comparing outcome between\n'...
    'patients with post-ictal change and those without is:\n'...
    'p = %1.2e\n'],p_post_outcome);

% age with post-ictal change
[p_post_age,info_post_age] = correlateClinically(postIcChange,age_all<=13,'bin','bin',0);
fprintf(['The p-value for chi squared comparing age between\n'...
    'patients with post-ictal change and those without is:\n'...
    'p = %1.2e\n'],p_post_age);
    
% correlate post-ictal change with temporal lobe onset
[p_post_lobe,info_post_lobe] = correlateClinically(postIcChange,temp_lobe_all,'bin','bin',0);
fprintf(['The p-value for chi squared comparing temporal vs non-temporal lobe between\n'...
    'patients with post-ictal change and those without is:\n'...
    'p = %1.2e\n'],p_post_lobe);


%% Table of all stuff

p_hour_all_t = getPText(p_hour_all);
p_pre_all_t = getPText(p_pre_all);
p_post_all_t = getPText(p_post_all);

T = table(names',p_hour_all_t,p_pre_all_t,p_post_all_t);


%% Table of hour long changes
table(char(names'),num2str(hour_chi2_all','%1.1f\n'),num2str(hour_dof_all','%d\n'),char(p_hour_all_t),'VariableNames',{'Name','Chi2','dof','P'})

%% Table of pre-ictal changes
table(char(names'),num2str(pre_chi2_all','%1.1f\n'),num2str(pre_dof_all','%d\n'),char(p_pre_all_t),'VariableNames',{'Name','Chi2','dof','P'})

%% Table of post-ictal changes
table(char(names'),num2str(post_chi2_all','%1.1f\n'),num2str(post_dof_all','%d\n'),char(p_post_all_t),'VariableNames',{'Name','Chi2','dof','P'})



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
    
    %{
    [p,tbl_kw,stats_kw] = kruskalwallis([SL_post_all,SL_inter_all,...
        SL_pre_all],[],'off');
    fprintf(['The p-value for change in SL comparing'...
        'pre/post/inter by K-W is:\n%1.2e\n'],p);
        %}
        
    [p,tbl,stats] = friedman([SL_post_all,SL_inter_all,...
        SL_pre_all],1,'off');
fprintf(['The Friedman test comparing SL is:\n'...
    'p = %1.2e, Q = %1.1f\n'],p,tbl{2,5});
    
    
    sig_post_move = find(p_post_all < 0.05/length(whichPts));
    names_post_move = names(sig_post_move)';
    SL_post_move = SL_post_all(sig_post_move);
    SL_other_move = SL_other_all(sig_post_move);
    SL_p_move = SL_p_post_all(sig_post_move);
    
    table(names_post_move,SL_post_move,SL_other_move,SL_p_move)
    
end

end



function p_text_cell = getPText(p_array)

    p_text_cell = cell(length(p_array),1);
    for i = 1:length(p_array)
        if p_array(i) < 0.0009
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

