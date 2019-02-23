function [n_grids_all,n_strips_all,n_depths_all] = getElectrodeInfo(pt,cluster,whichPts,printStuff)


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
end

n_strips_all = [];
n_depths_all = [];
n_grids_all = [];
n_elecs_all = [];

for whichPt = whichPts
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
    
    n_depths_all = [n_depths_all;n_depths];
    n_strips_all = [n_strips_all;n_strips];
    n_grids_all = [n_grids_all;n_grids];
    n_elecs_all = [n_elecs_all;n_elecs];
    
end

if printStuff == 1

fprintf('There is an average of %1.1f electrodes (range %d-%d).\n',...
    mean(n_elecs_all),min(n_elecs_all),max(n_elecs_all));

fprintf('There is an average of %1.1f grids (range %d-%d).\n',...
    mean(n_grids_all),min(n_grids_all),max(n_grids_all));

fprintf('There is an average of %1.1f strips (range %d-%d).\n',...
    mean(n_strips_all),min(n_strips_all),max(n_strips_all));

fprintf('There is an average of %1.1f depths (range %d-%d).\n',...
    mean(n_depths_all),min(n_depths_all),max(n_depths_all));

end

end