function cluster = getClusters(pt,whichPts)

%{
This function clusters spikes using k-means clustering based on their XYZ
coordinates of their spatial location. It also outputs sample spike
sequences for each cluster.
%}

%% Parameters
% Should we save the cluster struct?
saveStruct = 0;

% Merge with existing cluster struct (if 0, it will overwrite it)?
merge = 0; 

% get optimal cluster numbers? Set this to 1 if, instead of running the
% clustering algorithm, you instead want to run the algorithm to detect the
% optimal cluster number.
clustOpt = 1; 

% Plot example sequences? Should be 1 if doing as part of regular pipeline
doLongPlots = 1;

removeTies = 1; % I DO remove ties for every part of analysis of the spike paper
% My primary rationale for removing ties is that when I tested this for a
% single patient (HUP106), when I don't remove ties there is a third
% cluster that appears that is basically all artifact. This goes
% away when I remove ties.

% Get all spikes (not just lead spikes). This should be 1.
allSpikes = 1;

if isempty(whichPts) == 1
    for i = 1:length(pt)
        if isempty(pt(i).seq_matrix) == 0
            whichPts = [whichPts,i];
        end
    end
elseif whichPts == 100
    whichPts = [1,4,6,7,8,9,12,14,15,16,17,18,19,20,22,24,25,27,30,31];
end


%% Optimal cluster numbers
% These were chosen by elbow plot
if removeTies == 0
    n_clusters(1) = 3; % HUP64
    n_clusters(2) = 2; %HUP065
    n_clusters(3) = 4; %HUP68
    n_clusters(4) = 2; %HUP70
    n_clusters(5) = 2; %HUP73
    n_clusters(6) = 2; %HUP74
    n_clusters(7) = 2; %HUP75
    n_clusters(8) = 2; %HUP78
    n_clusters(9) = 3; %HUP080
    n_clusters(10) = 3; %HUP82
    n_clusters(11) = 3; %HUP83
    n_clusters(12) = 2; %HUP86 
    n_clusters(13) = 3; %HUP87
    n_clusters(14) = 2; %HUP88 
    n_clusters(15) = 3; %HUP94 
    n_clusters(16) = 2; %HUP105
    n_clusters(17) = 3; %HUP106 
    n_clusters(18) = 2; %HUP107
    n_clusters(19) = 2; %HUP111A
    n_clusters(20) = 4; %HUP116
    n_clusters(21) = 2; %Study012
    n_clusters(22) = 3; %Study16 
    n_clusters(23) = 2; %Study17
    n_clusters(24) = 3; %Study19
    n_clusters(25) = 4; %Study20
    n_clusters(26) = 0; %study21 - skip, no electrodes
    n_clusters(27) = 4; %Study22
    n_clusters(28) = 0; %study23 - skip, no electrodes
    n_clusters(29) = 0; %study26 - skip, no electrodes
    n_clusters(30) = 3; %Study28
    n_clusters(31) = 3; %Study29
elseif removeTies == 1
    n_clusters(1) = 3; % HUP64
    n_clusters(2) = 2; %HUP065
    n_clusters(3) = 4; %HUP68
    n_clusters(4) = 3; %HUP70
    n_clusters(5) = 4; %HUP73
    n_clusters(6) = 2; %HUP74
    n_clusters(7) = 2; %HUP75
    n_clusters(8) = 2; %HUP78
    n_clusters(9) = 3; %HUP080
    n_clusters(10) = 3; %HUP82
    n_clusters(11) = 3; %HUP83
    n_clusters(12) = 3; %HUP86
    n_clusters(13) = 3; %HUP87
    n_clusters(14) = 4; %HUP88
    n_clusters(15) = 4; %HUP94
    n_clusters(16) = 2; %HUP105
    n_clusters(17) = 2; %HUP106
    n_clusters(18) = 2; %HUP107
    n_clusters(19) = 2; %HUP111A
    n_clusters(20) = 4; %HUP116
    n_clusters(21) = 4; %Study012
    n_clusters(22) = 4; %Study16
    n_clusters(23) = 2; %Study17
    n_clusters(24) = 3; %Study19
    n_clusters(25) = 4; %Study20
    n_clusters(26) = 0; %study21 - skip, no electrodes
    n_clusters(27) = 4; %Study22
    n_clusters(28) = 0; %study23 - skip, no electrodes
    n_clusters(29) = 0; %study26 - skip, no electrodes
    n_clusters(30) = 3; %Study28
    n_clusters(31) = 3; %Study29
end

% Save file location
[~,~,scriptFolder,resultsFolder,~] = fileLocations;
p1 = genpath(scriptFolder);
addpath(p1);
destFolder = [resultsFolder,'clustering/validation/'];
mkdir(destFolder)

if merge == 1 && exist([destFolder,'cluster.mat'],'file') ~= 0
    temp = load([destFolder,'cluster.mat']);
    cluster = temp.cluster;
end

names = {};
all_sse = [];

for whichPt = whichPts
    
    % These are not actual patients
    if whichPt == 26 || whichPt == 28 || whichPt == 29
        continue
    end
    
    names = [names;pt(whichPt).name];
    
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
    
   
    %% Get all spike times and channels
    
    % all spike channels
    all_spikes = [];
    
    % all spike times
    all_times_all = [];
    
    % which sequence they're in
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
    
    %{
    I did not adjust seq_matrix after the ictal spike removal step. This
    should only be an issue in the CInfluence script, where I address it
    separately. Some of the sequences will not be accessed because I
    removed them from the sequence index. For example, if the indices
    1000-1010 are all ictal, and those belong to sequences 500 and 501,
    then those sequence indices will go away, and so I will not access
    those sequences for my example sequences.
    %}
    
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
        
        if 0
        figure
        plot(1:10,SSE)
        pause
        close(gcf)
        end
        
        if 0
            % Silhouette method
            all_ES = [];
            for i = 1:10
                tic
                E_S = evalclusters(all_locs,'kmeans','silhouette','klist',[1:10]);
                t = toc;
                all_ES = [all_ES;E_S.OptimalK,t];
            end

            % Gap method
            all_EG = [];
            for i = 1:10
                tic
                E_G = evalclusters(all_locs,'kmeans','gap','KList',[1:10]);
                t = toc;
                all_EG = [all_EG;E_G.OptimalK,t];
            end

            all_ES
            all_EG

            error('look\n');
        
        end
            
        all_sse = [all_sse;SSE'];
        
        
        
        
        
        continue
 
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
    
    if saveStruct == 1
        save([destFolder,'cluster.mat'],'cluster'); 
    end
        
    %% Get random group of sequences from each cluster and plot them
    if doLongPlots == 1
        rep_seq = [];
        for i = 1:n_clusters(whichPt)
            
            if allSpikes == 1
                % the indices of the spikes in this cluster
                spike_clust = find(idx==i); 

                % Take 50 of these indices randomly
                whichSpikes = spike_clust(randperm(length(spike_clust),...
                    min(50,length(spike_clust))));

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
    
  
end



if clustOpt == 1
    %% Plot the elbow plots for all patients
    n_clusters=n_clusters(whichPts);


    figure
    set(gcf,'position',[821 0 1113 800]);
    [ha, pos] = tight_subplot(4, 5, [0.05 0.03], [0.08 0.04], [0.04 0.01])
    for i = 1:length(names)
        axes(ha(i))
        plot(1:10,all_sse(i,:),'k','linewidth',2)
        hold on
        if i == 1
            plot(n_clusters(i),all_sse(i,n_clusters(i)),'marker','s',...
                'MarkerFaceColor','r','markeredgecolor','r','markersize',10)
            text(n_clusters(i) - 0.8, all_sse(i,n_clusters(i)) - 3.5e4,...
                sprintf('%d',n_clusters(i)),'fontsize',20);
        else
            plot(n_clusters(i),all_sse(i,n_clusters(i)),'marker','s',...
                'MarkerFaceColor','r','markeredgecolor','r','markersize',10)
            text(n_clusters(i) - 0.8, all_sse(i,n_clusters(i)) - 7e4,...
                sprintf('%d',n_clusters(i)),'fontsize',20);
        end
        title(sprintf('%s',names{i}))
        xlim([1,10])
        yticklabels([])

        xticks([1 10])


        if i == 18
            xlabel('Cluster number')
        end

        if i==11
            ylabel('                              Sum squared error')
        end
        set(gca,'fontsize',20)
    end
    print([destFolder,'elbows'],'-depsc')


end

end