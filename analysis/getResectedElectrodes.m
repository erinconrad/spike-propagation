function elecs = getResectedElectrodes(pt,whichPts)


% Save file location
[~,~,scriptFolder,resultsFolder,~] = fileLocations;
destFolder = [resultsFolder,'clustering/stats/'];
p1 = genpath(scriptFolder);
addpath(p1);
mkdir(destFolder);

%% Define which patients I am doing

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




end