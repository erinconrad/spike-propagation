function makeClinicalTables(pt,cluster,whichPts,doAll)
good = [];
if isempty(whichPts) == 1
    for i = 1:length(pt)
        if isempty(pt(i).seq_matrix) == 0
            if doAll == 1
                whichPts = [whichPts,i];
                if length(cluster(i).bad_cluster) < cluster(i).k
                    good = [good;1];
                else
                    good = [good;0];
                end
            else
                if length(cluster(i).bad_cluster) < cluster(i).k
                    whichPts = [whichPts,i];
                end
            end
        end
    end
    
    % I should be doing all 20 of these patients
    if isequal(whichPts,[1,4,6,7,8,9,12,14,15,16,17,18,19,20,22,24,25,27,30,31]) == 0
      %  error('Warning, not doing correct patients!\n');
    end
end

name = cell(length(whichPts),1);
ageOnset = cell(length(whichPts),1);
ageSurg = cell(length(whichPts),1);
sex = cell(length(whichPts),1);
soz = cell(length(whichPts),1);
path = cell(length(whichPts),1);
outcome = cell(length(whichPts),1);
grids = cell(length(whichPts),1);
depths = cell(length(whichPts),1);
strips =cell(length(whichPts),1);

count = 0;
for whichPt = whichPts
    count = count + 1;
    name{count} = pt(whichPt).name;
    ageOnset{count} = pt(whichPt).clinical.ageOnset;
    ageSurg{count} = pt(whichPt).clinical.ageSurgery;
    sex{count} = pt(whichPt).clinical.sex;
    soz{count} = pt(whichPt).clinical.seizureOnset;
    path{count} = pt(whichPt).clinical.pathology;
    %outcome{count} = getILAE(getOutcome(pt,whichPt));
    
    outcome{count} = pt(whichPt).clinical.outcome;
    
    n_elecs = length(pt(whichPt).channels);
    n_strips = 0;
    n_depths = 0;
    n_grids = 0;
    for i = 1:length(pt(whichPt).electrodeData.electrodes)
        if strcmp(pt(whichPt).electrodeData.electrodes(i).type,'S') == 1
            n_strips = n_strips + 1;
        elseif strcmp(pt(whichPt).electrodeData.electrodes(i).type,'G') == 1
            n_grids = n_grids + 1;
        elseif strcmp(pt(whichPt).electrodeData.electrodes(i).type,'D') == 1
            n_depths = n_depths + 1;
        else
            error('What\n');
        end
        
    end
    
    grids{count} = num2str(n_grids);
    strips{count} = num2str(n_strips);
    depths{count} = num2str(n_depths);
end

T = table(char(name),char(sex),char(ageOnset),char(ageSurg),char(soz),...
    (path),char(outcome),char(grids),char(strips),char(depths))

%% Compare good and bad for each column of the table

% Sex - bin/bin - chi squared
males = strcmp(sex,'M');
[tbl,chi2,p,labels] = crosstab(males,good)
fprintf('Chi2 test for men and women vs good and bad: p = %1.3f\n',p);

% Outcome - number/binary - wilcoxon rank sum
outcome_num = zeros(length(outcome),1);
for i = 1:length(outcome)
    t = outcome{i};
    outcome_num(i) = str2num(t(end));
end
[p,h,stats] = ranksum(outcome_num(good==1),outcome_num(good==0));
fprintf('Wilcoxon rank sum for outcomes for good detector vs bad: p = %1.3f\n',p);
[p, h, w_out] = ranksum_erin(outcome_num(good==1),outcome_num(good==0));
fprintf('Wilcoxon rank sum U = %1.1f\n\n',w_out);
test_stat = getStandardStats(outcome_num(good==1),outcome_num(good==0),'rs')
fprintf('Wilcoxon rank sum U = %1.1f\n\n',test_stat);

% Age at surg - num/bin - Wilcoxon rank sum
age_num = zeros(length(ageSurg),1);
for i = 1:length(ageSurg)
    t = ageSurg{i};
    if contains(t,'-') == 1
        a = regexp(t,'-');
        num1 = str2num(t(1:a-1));
        num2 = str2num(t(a+1:end));
        age_num(i) = mean([num1,num2]);
    elseif contains(t,'?') == 1
        age_num(i) = nan;
    else
        age_num(i) = str2num(t);
    end
end
[p,h,stats] = ranksum(age_num(good==1&isnan(age_num)==0),age_num(good==0&isnan(age_num)==0));
fprintf('Wilcoxon rank sum for age for good detector vs bad: p = %1.3f\n',p);

% Number of electrodes - num/bin - Wilcoxon rank sum
num_elecs = zeros(length(grids),1);
for i = 1:length(grids)
    grid_t = str2num(grids{i});
    strip_t = str2num(strips{i});
    depth_t = str2num(depths{i});
    num_elecs(i) = grid_t + strip_t + depth_t;
end
[p,h,stats] = ranksum(num_elecs(good==1),num_elecs(good==0));
fprintf('Wilcoxon rank sum for num elecs for good detector vs bad: p = %1.3f\n',p);

% TL vs not TL - bin/bin - chi2
tl = zeros(length(soz),1);
for i = 1:length(tl)
    if contains(soz{i},'T') == 1
        tl(i) = 1;
    else
        tl(i) = 0;
    end
end
[tbl,chi2,p,labels] = crosstab(tl,good)
fprintf('Chi2 test for TL and non-TL vs good and bad: p = %1.3f\n',p);


%% Get sex info
fprintf('There are %d males.\n',sum(strcmp(sex,'M')));

%% Get age info
a = ageSurg;
a{20} = nan;
a{18} = 25;
for i = 1:length(a)
    if isa(a{i},'char') == 1
        a{i} = str2double(a{i});
    end
end

b = cell2mat(a);
fprintf('The mean age is %1.1f, range %d-%d. %d peds.\n',...
    nanmean(b),min(b),max(b),sum(b<18));

end