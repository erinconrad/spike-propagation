function getClusters(pt,whichPts)

% NEED TO FIGURE OUT WHY THIS IS DIFFERENT FROM ORIGINAL

% Parameters
clustOpt = 1; % get optimal cluster numbers
doPlots = 1;
doLongPlots = 1;
removeTies = 1;

% Optimal cluster numbers
n_clusters(3) = 4; %HUP68, 2 by silhouette
n_clusters(4) = 2; %HUP70
n_clusters(8) = 3; %HUP78
n_clusters(9) = 3; %HUP080
n_clusters(12) = 3; %HUP86
n_clusters(17) = 2; %HUP106
n_clusters(18) = 4; %HUP107
n_clusters(19) = 3; %HUP111A
n_clusters(20) = 3; %HUP116
n_clusters(22) = 4; %Study16
n_clusters(24) = 3; %Study19
n_clusters(25) = 4; %Study20
n_clusters(27) = 3; %Study22
n_clusters(30) = 4; %Study28
n_clusters(31) = 3; %Study29

% Save file location
[~,~,~,resultsFolder,~] = fileLocations;
destFolder = [resultsFolder,'clustering/validation/'];
mkdir(destFolder)

for whichPt = whichPts
    
    % Patient parameters
    fprintf('Doing %s\n',pt(whichPt).name);
    locs = pt(whichPt).electrodeData.locs;
    szTimes = pt(whichPt).newSzTimes;
    saveFolder = [destFolder,pt(whichPt).name,'/'];
    mkdir(saveFolder);
    
    %% Get spike times and locations
    
    % Get all sequences
    seq_matrix = pt(whichPt).seq_matrix;
    
    % Get all spike times and channels
    all_spikes = [];
    all_times_all = [];
    seq_index = [];
    for i = 1:size(seq_matrix,2)
        nonan = find(~isnan(seq_matrix(:,i)));
        times = seq_matrix(nonan,i);
        
        % Resort by time
        [~,I] = sort(times);
        nonan = nonan(I);
        
        all_spikes = [all_spikes;nonan];
        all_times_all = [all_times_all;seq_matrix(nonan,i)];
        
        
        seq_index = [seq_index;i*ones(length(nonan),1)];
    end
    
    % Remove ictal spikes and spikes within 1 minute of seizure start
    t = find(any(all_times_all >= (szTimes(:,1)-repmat(60,size(szTimes,1),1))' ...
        & all_times_all <= szTimes(:,2)',2));
    all_times_all(t) = [];
    all_spikes(t) = [];
    fprintf('Removed %d ictal spikes \n',length(t));
    
    % Get locations of spikes
    all_locs = locs(all_spikes,:);
    
    % Fill cluster structure
    pt(whichPt).cluster.all_spikes = all_spikes;
    pt(whichPt).cluster.all_locs = all_locs;
    pt(whichPt).cluster.all_times_all = all_times_all;
    
    %% Determine optimal number of clusters
    if clustOpt == 1
        
        % Elbow method
        SSE = zeros(10,1);
        for k = 1:10
            
            SSE_temp = zeros(30,1);
            
            % Do clustering algorithm 30 times, take the best
            for j = 1:30
                % Cluster
                [idx_test,C_test] = kmeans(all_locs,k);
                
                 % Get SSE
                 for i = 1:k
                    % sum of square differences between each observation
                    % and its cluster's centroid
                    % CHECK ME
                    SSE_temp(j) = SSE_temp(j) + ...
                        sum(sum((all_locs(idx_test == i,:) - ...
                        repmat(C_test(i,:),...
                        size(all_locs(idx_test == i,:),1),1)).^2));
                 end
            end  
            SSE(k) = min(SSE_temp);
        end
        figure
        plot(1:10,SSE)
        
        % Silhouette method
        %E_S = evalclusters(all_locs,'kmeans','silhouette','klist',[1:10]);
        
        % Gap method
        %E_G = evalclusters(all_locs,'kmeans','gap','KList',[1:10]);
 
    end
    
    
end


end