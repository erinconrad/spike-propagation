function cluster = getClusters(pt,whichPts)

%% Parameters
% Should we skip patients that are already done and merge the new patients
% with the existing cluster?
merge = 1; 
allSpikes = 1;
clustOpt = 0; % get optimal cluster numbers
doPlots = 0;
doLongPlots = 1;
removeTies = 0;


%% Optimal cluster numbers
n_clusters(3) = 4; %4; %HUP68, 2 by silhouette
n_clusters(4) = 3; %2; %HUP70
n_clusters(8) = 2; %3; %HUP78
n_clusters(9) = 6; %3; %HUP080
n_clusters(12) = 3; %3; %HUP86
n_clusters(17) = 2; %2; %HUP106
n_clusters(18) = 2; %4; %HUP107
n_clusters(19) = 2; %3; %HUP111A
n_clusters(20) = 4; %3; %HUP116
n_clusters(22) = 4; %4; %Study16
n_clusters(24) = 3; %3; %Study19
n_clusters(25) = 4; %4; %Study20
n_clusters(27) = 4; %3; %Study22
n_clusters(30) = 3; %4; %Study28
n_clusters(31) = 3; %3; %Study29

% Save file location
[~,~,scriptFolder,resultsFolder,~] = fileLocations;
p1 = genpath(scriptFolder);
addpath(p1);
destFolder = [resultsFolder,'clustering/validation/'];
mkdir(destFolder)

if merge == 1
    temp = load([destFolder,'cluster.mat']);
    cluster = temp.cluster;
end

for whichPt = whichPts
    
    
    
    
    cluster(whichPt).name = pt(whichPt).name;
    cluster(whichPt).allSpikes = allSpikes;
    
    
    % Patient parameters
    fprintf('Doing %s\n',pt(whichPt).name);
    locs = pt(whichPt).electrodeData.locs(:,2:4);
    szTimes = pt(whichPt).newSzTimes;
    saveFolder = [destFolder,pt(whichPt).name,'/'];
    
    
    % I can skip doing a patient if I already did them
    if merge == 1 && exist(saveFolder,'dir') == 7
        skipFlag = 0;
        for i = 1:length(cluster)
            if strcmp(pt(whichPt).name,cluster(i).name) == 1
                fprintf('Already did %s, skipping\n',...
                    cluster(i).name);
                skipFlag = 1;
                break
            end
        end
        if skipFlag == 1
            continue
        end
    end
    
    mkdir(saveFolder);
    
    %% Get spike times and locations
    
    % Get all sequences
    seq_matrix = pt(whichPt).seq_matrix;
    
    % Remove sequences with too many ties
    if removeTies == 1
        keep = ones(size(seq_matrix,2),1);
        for s = 1:size(seq_matrix,2)
           curr_seq = seq_matrix(:,s);
           nonans = curr_seq(~isnan(curr_seq));
           norepeats = unique(nonans);
           if length(norepeats) < 0.5*length(nonans)
               keep(s) = 0;
           end
        end
        seq_matrix(:,keep==0) = [];
        fprintf(['%s had %d sequences (%1.2f of all sequences) deleted'...
        'for having >50 percent ties\n%d sequences remain\n'],...
        pt(whichPt).name,sum(keep == 0),sum(keep == 0)/length(keep),sum(keep==1));

    end
    
    % Get first channel times and locations
    [lead_times,lead] = min(seq_matrix,[],1);
    lead_locs = locs(lead,:);
    cluster(whichPt).lead_locs = lead_locs;
    cluster(whichPt).lead_times = lead_times;
    
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
    seq_index(t) = [];
    fprintf('Removed %d ictal spikes \n',length(t));
    
    % Get locations of spikes
    all_locs = locs(all_spikes,:);
    
    % Fill cluster structure
    cluster(whichPt).all_spikes = all_spikes;
    cluster(whichPt).all_locs = all_locs;
    cluster(whichPt).all_times_all = all_times_all;
    cluster(whichPt).seq_index = seq_index;
    
    %% Determine optimal number of clusters
    if clustOpt == 1
        
        % Elbow method
        SSE = zeros(10,1);
        for k = 1:10
            
            SSE_temp = zeros(30,1);
            
            % Do clustering algorithm 30 times, take the best
            for j = 1:30
                % Cluster
                if allSpikes == 1
                    [idx_test,C_test] = kmeans(all_locs,k);
                else
                    [idx_test,C_test] = kmeans(lead_locs,k);
                end
                
                 % Get SSE
                 for i = 1:k
                    % sum of square differences between each observation
                    % and its cluster's centroid
                    % CHECK ME
                    if allSpikes == 1
                        SSE_temp(j) = SSE_temp(j) + ...
                            sum(sum((all_locs(idx_test == i,:) - ...
                            repmat(C_test(i,:),...
                            size(all_locs(idx_test == i,:),1),1)).^2));
                    else
                        SSE_temp(j) = SSE_temp(j) + ...
                            sum(sum((lead_locs(idx_test == i,:) - ...
                            repmat(C_test(i,:),...
                            size(lead_locs(idx_test == i,:),1),1)).^2));
                    end
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
    
    %% Do clustering algorithm
    for j = 1:30 % do it 30 times and take the best result
        if allSpikes == 1
            [idx_all{j},C_all{j},sumd_all{j},D_all{j}] = ...
            kmeans([all_locs],n_clusters(whichPt));
            metric(j) = sum(sumd_all{j});
        else
            [idx_all{j},C_all{j},sumd_all{j},D_all{j}] = ...
            kmeans([lead_locs],n_clusters(whichPt));
            metric(j) = sum(sumd_all{j});
        end
        
        %{
        (for testing purposes. This gives the same result as the above
        metric).
        alt_metric(j) = 0;
        for i = 1:n_clusters(whichPt)
            alt_metric(j) = alt_metric(j) + ...
                sum(sum((all_locs(idx_all{j} == i,:) - ...
                repmat(C_all{j}(i,:),...
                size(all_locs(idx_all{j} == i,:),1),1)).^2));
        end
        %}
    end
    
    % Take the result of the clustering algorithm that worked best
    [~,minidx] = min(metric);
    idx = idx_all{minidx};
    C = C_all{minidx};
    D = D_all{minidx};
    
    cluster(whichPt).idx = idx;
    cluster(whichPt).C = C;
    cluster(whichPt).D = D;
    cluster(whichPt).k = n_clusters(whichPt);
    cluster(whichPt).bad_cluster = [];
    
    save([destFolder,'cluster.mat'],'cluster'); 
    
    %% Get random group of sequences from each cluster and plot them
    if doLongPlots == 1
        rep_seq = [];
        for i = 1:n_clusters(whichPt)
            
            if allSpikes == 1
                % the indices of the spikes in this cluster
                spike_clust = find(idx==i); 

                % Take 50 of these indices randomly
                whichSpikes = spike_clust(randperm(length(spike_clust),50));

                % Get the sequences these spikes belong to
                whichSeqs = seq_index(whichSpikes);
                rep_seq{i} = seq_matrix(:,whichSeqs);
            else
                [sortedD,I] = sort(D(:,i));
                spike_clust = find(idx==i);
                rep_seq{i} = seq_matrix(:,I(randperm(length(I),12)));
            end
            
            % Plot the sequences
            outputFolder = [saveFolder,sprintf('seqs_cluster_%d',i),'/'];
            showSequences(pt,whichPt,rep_seq{i},[],0,outputFolder)
            
            % Plot gifs
            info(i).cluster = i;
            info(i).name = pt(whichPt).name;
            info(i).outputFile = [outputFolder,'cluster_',sprintf('%d',i),'.gif'];
            movieSeqs(rep_seq{i}(:,1:10),locs,C(i,:),info(i));
        end
    end
    
    %% Plot the changing cluster over time
    if doPlots == 1
        % Assign the sequence a color based on its cluster index
        colors = [0 0 1;1 0 0;0 1 0; 0.5 0.5 1; 1 0.5 0.5; 0.5 1 0.5; 0.4 0.7 0.4];
        c_idx = zeros(size(idx,1),3);
        for i = 1:length(idx)
           c_idx(i,:) = colors(idx(i),:); 
        end
        
        figure
        set(gcf,'Position',[50 100 1200 1200])
        
        %% Subplot 1: Plot the x, y, z over time
        subplot(3,1,1)
        ttext = {'x','y','z'};
        toAdd = 0;
        if allSpikes == 1
            plot_thing = all_locs;
            plot_times = all_times_all;
        else
            plot_thing = lead_locs;
            plot_times = lead_times';
        end
        
        for i = 1:3
            % Plot each coordinate over time, offset from each other
            scatter(plot_times/3600,plot_thing(:,i)+...
                repmat(toAdd,size(plot_thing,1),1),20,c_idx)
            hold on
            text(plot_times(1)/3600-0.3,toAdd+median(plot_thing(:,i)),...
                sprintf('%s',ttext{i}),'FontSize',30);
            if i ~=3
                % Define the offset for each coordinate
                toAdd = toAdd + 10+(max(plot_thing(:,i)) - ...
                    min(plot_thing(:,i+1)));
            end
        end
        
        % Plot the seizure times
        for j = 1:size(szTimes,1) 
            yl = ylim;
            plot([szTimes(j,1) szTimes(j,1)]/3600,yl,'k','LineWidth',2);
        end
        
        set(gca,'ytick',[]);
        xlim([plot_times(1)/3600-1 plot_times(end)/3600+1])
        title(sprintf('X, y, z coordinates of all spikes for %s',...
            pt(whichPt).name));
        set(gca,'FontSize',15);
        
        %% Subplot 2: Plot the cluster distribution over time
        
        % Get the times for spikes in each cluster
        for i = 1:n_clusters(whichPt)
            clust{i} = plot_times(idx == i);
        end
        
        
        window = 3600;
        
        % NEED TO CHECK THIS SCRIPT - ALSO IT's SUUUPER SLOW
        %{
        [sum_c,sum_times] = movingSumCounts(clust,plot_times,window);
        totalSum = zeros(1,size(sum_times,2));
        for i = 1:n_clusters(whichPt)
            totalSum = totalSum + sum_c(i,:);
        end
        prop_c = sum_c./totalSum;
        %}
        
        
        % alternate method - try non-moving window
        
        % enough bins so that it's essentially 10 minute windows
        nbins = round((max(plot_times)-min(plot_times))/(window/3));
        [Y,E] = discretize(plot_times,nbins);
        new_times = E(2:end);
        new_counts = zeros(nbins,n_clusters(whichPt));
        for bb = 1:nbins
            for k = 1:n_clusters(whichPt)
                new_counts(bb,k) = sum(Y==bb & idx==k);
            end
        end
        new_prop = new_counts./sum(new_counts,2);
        sum_times = new_times;
        prop_c = new_prop';
        
        subplot(3,1,2)
        pl = zeros(n_clusters(whichPt),1);
        for i = 1:n_clusters(whichPt)
            pl(i)= plot(sum_times/3600,prop_c(i,:),...
                'color',colors(i,:),'LineWidth',2);
        hold on
        end
        
        for j = 1:size(szTimes,1) 
            yl = ylim;
            plot([szTimes(j,1) szTimes(j,1)]/3600,yl,'k','LineWidth',2);
        end
        xlim([plot_times(1)/3600-1 plot_times(end)/3600+1])
        title(sprintf(['Proportion of sequences in given cluster, moving'...
        ' average %d s, %s'],window,pt(whichPt).name));
        set(gca,'FontSize',15);
        
        %% Subplot 3: Plot locations of centroids
        subplot(3,1,3)
        scatter3(locs(:,1),locs(:,2),locs(:,3),60,'k');
        hold on
        for k = 1:size(C,1)
            scatter3(C(k,1),C(k,2),C(k,3),60,colors(k,:),'filled');
        end
        title(sprintf('Spike location centroids for each cluster for %s',...
            pt(whichPt).name))
        set(gca,'FontSize',15);
        set(gca,'xticklabel',[])
        set(gca,'yticklabel',[])
        set(gca,'zticklabel',[])
        
    end
    
end


end