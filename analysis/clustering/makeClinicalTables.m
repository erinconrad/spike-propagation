function makeClinicalTables(pt,cluster,whichPts)

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
    
    grids{count} = n_grids;
    strips{count} = n_strips;
    depths{count} = n_depths;
end

T = table(name,sex,ageOnset,ageSurg,soz,path,outcome,grids,strips,depths);


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