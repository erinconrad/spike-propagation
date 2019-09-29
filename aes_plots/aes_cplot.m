function aes_cplot(pt,cluster,whichPt)

% Show absolute number?
show_abs = 0;
show_sz = 0;

% Save file location
[~,~,scriptFolder,resultsFolder,~] = fileLocations;
destFolder = [resultsFolder,'clustering/plots/'];
mkdir(destFolder);
addpath(genpath(scriptFolder))

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

all_times_all = cluster(whichPt).all_times_all; % all spike times
all_spikes = cluster(whichPt).all_spikes; % all spike channels
all_locs = cluster(whichPt).all_locs;
k = cluster(whichPt).k; % the number of clusters
idx = cluster(whichPt).idx; % the cluster index for every spike
C = cluster(whichPt).C; % the centroids of the clusters
bad_cluster = cluster(whichPt).bad_cluster; % which clusters are bad

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
all_locs(bad_idx,:) = [];
idx(bad_idx) = [];
clusters = 1:k; clusters(bad_cluster) = [];
C(bad_cluster,:) = [];

colors = [0 0 1;1 0 0;0 1 0];
c_idx = zeros(size(idx,1),3);

plot_thing = all_locs;
plot_times = all_times_all;
% Get the times for spikes in each cluster
for i = 1:length(clusters)
    clust{i} = plot_times(idx == clusters(i));
end
window = 3600;
sparse_time = plot_times-min(plot_times);

[sum_c,sum_times] = movingSumCounts(clust,plot_times,window);
totalSum = zeros(1,size(sum_times,2));
for i = 1:length(clusters)
    totalSum = totalSum + sum_c(i,:);
end
prop_c = sum_c./totalSum;

for i = 1:length(clusters)  
   c_idx(idx==clusters(i),:) = ...
       repmat(colors(i,:),sum(idx==clusters(i)),1);
end

possibleText = {'Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5','Cluster 6'};
textLeg = possibleText(1:size(C,1));


if show_abs == 1
    figure
    set(gcf,'Position',[71 241 1237 538])
    [ha, ~] = tight_subplot(2, 1, [.09 .03],[.14 .07],[.01 .01]);

    axes(ha(1));
    pl = zeros(length(clusters),1);
    for i = 1:length(clusters)
        pl(i)= plot((sum_times-min(plot_times))/3600,sum_c(i,:),...
            'color',colors((i),:),'LineWidth',2);
    hold on
    end
    xlim([sparse_time(1)/3600-1 sparse_time(end)/3600+1])
    %ylabel('Number of spikes')
    yticklabels([])
    xticklabels([])
    title(sprintf(['Number of spikes in given cluster '...
    '(moving average)']))
    l = legend(pl,textLeg,'location','northwest');
    set(gca,'FontSize',25);

    axes(ha(2));
    pl = zeros(length(clusters),1);
    for i = 1:length(clusters)
        pl(i)= plot((sum_times-min(plot_times))/3600,prop_c(i,:),...
            'color',colors((i),:),'LineWidth',2);
    hold on
    end
    %ylabel('Proportion of spikes')
    yticklabels([])
    xlim([sparse_time(1)/3600-1 sparse_time(end)/3600+1])
    xlabel('Time (hours)');
    title(sprintf(['Proportion of spikes in given cluster '...
    '(moving average)']))
    set(gca,'FontSize',25);

else
    
    figure
    set(gcf,'Position',[71 241 1237 300])
    [ha, ~] = tight_subplot(1, 1, [.09 .03],[.24 .13],[.04 .01]);
    axes(ha(1))
    pl = zeros(length(clusters),1);
    for i = 1:length(clusters)
        pl(i)= plot((sum_times-min(plot_times))/3600,prop_c(i,:),...
            'color',colors((i),:),'LineWidth',2);
    hold on
    end
    %ylabel('Proportion of spikes')
    %yticklabels([])
    xlim([sparse_time(1)/3600-1 sparse_time(end)/3600+1])
    xlabel('Time (hours)');
    title(sprintf(['Proportion of spikes in given cluster '...
    '(moving average)']))
    set(gca,'FontSize',25);
    
    if show_sz == 1
        for j = 1:size(szTimes,1) 
            yl = ylim;
            plot([szTimes(j,1)-min(plot_times) szTimes(j,1)-min(plot_times)]/3600,...
                yl,'k','LineWidth',2);
        end
        
        if whichPt == 31
            annotation('textarrow',[0.22 0.285],[0.55 0.55],'String','Seizure',...
                'FontSize',25);

            annotation('textarrow',[0.5-.075 0.5],[0.55 0.55],'String','Seizure',...
                'FontSize',25);

            annotation('textarrow',[0.74-.075 0.744],[0.55 0.55],'String','Seizure',...
                'FontSize',25);
        end
    end
    
    
end


end