function CPlots(pt,cluster,whichPts)

%{ 

CPlots
This is my cleaned up file for getting plots on the cluster data 

%}    

% Save file location
[~,~,~,resultsFolder,~] = fileLocations;
destFolder = [resultsFolder,'clustering/plots/'];
mkdir(destFolder);

allCounts = [];
allPat = [];
allChunk = [];


for whichPt = whichPts
    
    fprintf('Doing %s\n',pt(whichPt).name);
    
    % Get patient parameters
    locs = pt(whichPt).electrodeData.locs(:,2:4);
    szTimes = pt(whichPt).newSzTimes;
    soz = pt(whichPt).newSOZChs;
    
    % Reorder seizure times if out of order
    oldSzTimes = szTimes;
    szTimes = sort(szTimes,1);
    if isequal(oldSzTimes,szTimes) == 0
        fprintf('WARNING!!! %s seizure times out of order\n',pt(whichPt).name);
    end
    
    % Combine nearly equal seizure times
    newIdx = 2;
    newSzTimes = [];
    newSzTimes(1,:) = szTimes(1,:);
    for j = 2:size(szTimes,1)
        if abs(szTimes(j,1)-szTimes(j-1,1)) < 10 && ...
                abs(szTimes(j,2)-szTimes(j-1,2))
           newIdx = newIdx - 1; 
           newSzTimes(newIdx,1) = min(szTimes(j,1),szTimes(j-1,1));
           newSzTimes(newIdx,2) = max(szTimes(j,2),szTimes(j-1,2));  
        else
           newSzTimes(newIdx,:) = szTimes(j,:);
        end
        newIdx = newIdx + 1;
    end
    
    if isequal(newSzTimes,szTimes) == 0
        fprintf('WARNING!!! %s had duplicate seizure times\n',pt(whichPt).name);
    end
    
    % Pull cluster info
    %{
    all_times_all = pt(whichPt).cluster.all_times_all; % all spike times
    all_spikes = pt(whichPt).cluster.all_spikes; % all spike channels
    all_locs = pt(whichPt).cluster.all_locs;
    k = pt(whichPt).cluster.k; % the number of clusters
    idx = pt(whichPt).cluster.idx; % the cluster index for every spike
    C = pt(whichPt).cluster.C; % the centroids of the clusters
    bad_cluster = pt(whichPt).cluster.bad_cluster; % which clusters are bad
    %}
    all_times_all = cluster(whichPt).all_times_all; % all spike times
    all_spikes = cluster(whichPt).all_spikes; % all spike channels
    all_locs = cluster(whichPt).all_locs;
    k = cluster(whichPt).k; % the number of clusters
    idx = cluster(whichPt).idx; % the cluster index for every spike
    C = cluster(whichPt).C; % the centroids of the clusters
    bad_cluster = cluster(whichPt).bad_cluster; % which clusters are bad
    allSpikes = cluster(whichPt).allSpikes; % am I taking all spikes???
    
    % Confirm that I do not have any ictal spikes
    t = find(any(all_times_all >= szTimes(:,1)' & all_times_all <= szTimes(:,2)',2));
    if isempty(t) == 0
        fprintf('WARNING: Remaining ictal spikes for %s!\n',pt(whichPt).name);
        all_times_all(t) = [];
        all_spikes(t) = [];
        all_locs(t,:) = [];
        idx(t) = [];
    end
    

    % Remove bad clusters
    bad_idx = find(ismember(idx,bad_cluster));
    all_times_all(bad_idx) = [];
    all_spikes(bad_idx) = [];
    all_locs(bad_idx,:) = [];
    idx(bad_idx) = [];
    clusters = 1:k; clusters(bad_cluster) = [];
    
    %% Do plots
    
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
    plot_thing = all_locs;
    plot_times = all_times_all;

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
    for i = 1:length(clusters)
        clust{i} = plot_times(idx == i);
    end


    window = 1800;

    % NEED TO CHECK THIS SCRIPT - ALSO IT's SUUUPER SLOW
    
    [sum_c,sum_times] = movingSumCounts(clust,plot_times,window);
    totalSum = zeros(1,size(sum_times,2));
    for i = 1:length(clusters)
        totalSum = totalSum + sum_c(i,:);
    end
    prop_c = sum_c./totalSum;
    %}


    % alternate method - try non-moving window

    % enough bins so that it's essentially 30 minute windows
    nbins = round((max(plot_times)-min(plot_times))/(window));
    [Y,E] = discretize(plot_times,nbins);
    new_times = E(2:end);
    new_counts = zeros(nbins,length(clusters));
    for bb = 1:nbins
        for k = 1:length(clusters)
            new_counts(bb,k) = sum(Y==bb & idx==k);
        end
    end
    new_prop = new_counts./sum(new_counts,2);
    %sum_times = new_times;
    %prop_c = new_prop';

    subplot(3,1,2)
    pl = zeros(length(clusters),1);
    for i = 1:length(clusters)
        pl(i)= plot(sum_times/3600,prop_c(i,:),...
            'color',colors(i,:),'LineWidth',2);
    hold on
    end

    for j = 1:size(szTimes,1) 
        yl = ylim;
        plot([szTimes(j,1) szTimes(j,1)]/3600,yl,'k','LineWidth',2);
    end
    xlim([plot_times(1)/3600-1 plot_times(end)/3600+1])
    title(sprintf(['Proportion of sequences in given cluster, '...
    '%d minute bins, %s'],window/60,pt(whichPt).name));
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