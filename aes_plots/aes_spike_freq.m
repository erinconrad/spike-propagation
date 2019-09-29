function aes_spike_freq(pt,cluster)

whichPt = 31;
circ_size = 1000;

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;

locs = pt(whichPt).electrodeData.locs(:,2:4);
chs = 1:size(locs,1);
szTimes = pt(whichPt).newSzTimes;
soz = pt(whichPt).newSOZChs;

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

%% Get cluster identity of each channel
for i = 1:3
    cl_locs{i} = unique(all_spikes(idx == i,:),'rows');
end


% Get spike counts in each channel
ch_counts = zeros(length(chs),1);
for k = 1:length(chs)
    ch_counts(k) = sum(all_spikes==k);
end

% Get maximum brightness
alpha_bright = linspace(1,0,max(ch_counts));

%% Plot
figure
set(gcf,'position',[10 10 900 900])
set(gcf,'color','white');
rc = scatter3(locs(:,1),locs(:,2),locs(:,3),circ_size,'r','filled');
%set(rc,'visible','off');
hold on
bc = scatter3(locs(:,1),locs(:,2),locs(:,3),circ_size,'b','filled');
%set(bc,'visible','off');
gc = scatter3(locs(:,1),locs(:,2),locs(:,3),circ_size,'g','filled');
%set(gc,'visible','off');
scatter3(locs(:,1),locs(:,2),locs(:,3),circ_size,'w','filled');
scatter3(locs(:,1),locs(:,2),locs(:,3),circ_size,'k','linewidth',2);

for k = 1:length(chs)
    
    % get the appropriate cluster
    for cl = 1:length(cl_locs)
        if ismember(k,cl_locs{cl}) == 1
            which_clust = cl;
        end
    end
    
    % Get appropriate color
    n_counts = ch_counts(k);
    
    if which_clust == 1
        this_col = [alpha_bright(max(1,n_counts)) ...
            alpha_bright(max(1,n_counts)) 1];
    elseif which_clust == 2
        this_col = [1 alpha_bright(max(1,n_counts)) alpha_bright(max(1,n_counts))];
    elseif which_clust == 3
        this_col = [alpha_bright(max(1,n_counts)) 1 ...
            alpha_bright(max(1,n_counts))];
    end
    
    scatter3(locs(k,1),locs(k,2),locs(k,3),circ_size,this_col,'filled');
    
    
end

if 1 == 1
    sc = scatter3(locs(soz,1),locs(soz,2),locs(soz,3),circ_size-50,'p',...
        'markerfacecolor','w','markeredgecolor',...
        'k','linewidth',2);

    l = legend([bc,rc,gc,sc],{'Cluster 1','Cluster 2','Cluster 3','Seizure onset zone'},...
        'location','northeast');
else
    l = legend([bc,rc,gc],{'Cluster 1','Cluster 2','Cluster 3'},...
        'location','northeast');
end
view(24,2.4)
pause(0.5);
for i = 1:length(l.EntryContainer.NodeChildren)
    l.EntryContainer.NodeChildren(i).Icon.Transform.Children.Children.Size = 20;
end
set(gca,'fontsize',25)
xticklabels([])
yticklabels([])
zticklabels([])
grid off

end