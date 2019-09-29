function aes_brain_gif(pt,cluster)

% This is for HUP078
whichPt = 8;
offset = [-3 27 7.7653];
window = 3600;

%% Get stuff
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile,other] = fileLocations;
addpath(other.gifti)
outputFolder = [resultsFolder,'pretty_plots/aes_gif/'];
mkdir(outputFolder)
other_file_out = [outputFolder,'elecs_aes.gif'];


locs = pt(whichPt).electrodeData.locs(:,2:4);

% Get transformation matrix to get new coordinate locations
A = makeNewElecData(pt,whichPt);
%offset = [-10 30 0]; % bs
locs = A*locs-offset;
soz = pt(whichPt).newSOZChs;
szTimes = pt(whichPt).newSzTimes;
chs = 1:size(locs,1);

%% Load gifti
brainFolder = '/Users/erinconrad/Desktop/residency stuff/R25/actual work/data/brains/';
giftiFolder = [brainFolder,pt(whichPt).name,'/'];
names = dir([giftiFolder,'*pial.gii']);
fname2 = names(1).name;
g = gifti([giftiFolder,fname2]);

%% Reorder seizure times if out of order
oldSzTimes = szTimes;
szTimes = sort(szTimes,1);
if isequal(oldSzTimes,szTimes) == 0
    fprintf('WARNING!!! %s seizure times out of order\n',pt(whichPt).name);
end

%% Combine nearly equal seizure times
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

%% Get cluster info
all_times_all = cluster(whichPt).all_times_all; % all spike times
all_spikes = cluster(whichPt).all_spikes; % all spike channels
all_locs = cluster(whichPt).all_locs;
k = cluster(whichPt).k; % the number of clusters
idx = cluster(whichPt).idx; % the cluster index for every spike
C = cluster(whichPt).C; % the centroids of the clusters
bad_cluster = cluster(whichPt).bad_cluster; % which clusters are bad

%% Confirm that I do not have any ictal spikes
t = find(any(all_times_all >= szTimes(:,1)' & all_times_all <= szTimes(:,2)',2));
if isempty(t) == 0
    fprintf('WARNING: Remaining ictal spikes for %s!\n',pt(whichPt).name);
    all_times_all(t) = [];
    all_spikes(t) = [];
    all_locs(t,:) = [];
    idx(t) = [];
end


%% Remove bad clusters
bad_idx = find(ismember(idx,bad_cluster));
all_times_all(bad_idx) = [];
all_locs(bad_idx,:) = [];
idx(bad_idx) = [];
clusters = 1:k; clusters(bad_cluster) = [];
C(bad_cluster,:) = [];
all_spikes(bad_idx) = [];
    
%% Get cluster identity of each channel
for i = 1:length(clusters)
    cl_locs{i} = unique(all_spikes(idx == clusters(i),:),'rows');
end


plot_times = all_times_all;

%% Find if data is broken into chunks (like for HUP080)

diff_times = diff(plot_times);
big_diff = diff_times > 1e5;

change_idx = find(big_diff);
group_idx = zeros(length(change_idx)+1,2);

if isempty(change_idx) == 1
    group_idx = [1 length(plot_times)];
else
    group_idx(1,:) = [1 change_idx(1)];
    for i = 2:size(group_idx,1)
        group_idx(i,1) = change_idx(i-1)+1;
        if i == size(group_idx,1)
            group_idx(i,2) = length(plot_times);
        else
            group_idx(i,2) = change_idx(i);
        end
    end
end

for i = 1:sum(big_diff) + 1
    
curr_idx = group_idx(i,:);
    
%% Define windows
nbins = ceil((plot_times(curr_idx(2))-plot_times(curr_idx(1)))/window);
window_edges = zeros(nbins,2);
for j = 1:nbins
    window_edges(j,1) = plot_times(curr_idx(1)) + window*(j-1);
    window_edges(j,2) = min(plot_times(curr_idx(1)) + window*(j),plot_times(curr_idx(2)));
end
E = [window_edges(:,1);window_edges(end,2)];
[Y] = discretize(plot_times(curr_idx(1):curr_idx(2)),E);
new_times = round(E(1:end)/3600);
new_times = new_times-min(new_times) + 1;
new_counts = zeros(nbins,length(clusters));

idx_to_look = idx(curr_idx(1):curr_idx(2));
for bb = 1:nbins
    for k = 1:length(clusters)
        new_counts(bb,k) = sum(Y==bb & idx_to_look==clusters(k));
    end
end
new_prop = new_counts./sum(new_counts,2);

%% Bin times - number in each channel over time
[Y] = discretize(plot_times(curr_idx(1):curr_idx(2)),E);
ch_counts = zeros(nbins,length(chs));
bin_times = zeros(nbins,2);
for bb = 1:nbins
    for k = 1:length(chs)
        ch_counts(bb,k) = sum(Y==bb & all_spikes(curr_idx(1):curr_idx(2))==k);
    end
    bin_times(bb,:) = [E(bb),E(bb+1)];
end


if 0
    %% Plot single time
    tt = 10;
    
    fig = figure;
    set(gcf,'color','white');
    
    p = plotGIFTI(g);
    hold on
    
    % Get counts per channel in this time
    alpha_lin_time = linspace(1,0,max(ch_counts(tt,:)));
    alpha_lin_time_yellow = [linspace(1,0.9290,max(ch_counts(tt,:))); ...
        linspace(1,0.540,max(ch_counts(tt,:))); linspace(1,0.1250,max(ch_counts(tt,:)))];
    
    scatter3(locs(:,1),locs(:,2),locs(:,3),350,'k','linewidth',2);
    hold on
    %view(75.6 - tt*3,-4.2)
    if whichPt == 31
        view(75.6,-4.2)
    elseif whichPt == 9
        view(104.4,14.2)
    elseif whichPt == 17
        view(-196.3,1.2)
    elseif whichPt == 8
        view(-120,-11);
    end
    
    
    for k = 1:length(chs)
        
        % get the appropriate cluster
        for cl = 1:length(cl_locs)
            if ismember(k,cl_locs{cl}) == 1
                which_clust = cl;
            end
        end
        
        % Get appropriate color
        n_counts = ch_counts(tt,k);
        
        %{
        if which_clust == 1
            this_col = [alpha_lin_time(max(1,n_counts)) ...
                alpha_lin_time(max(1,n_counts)) 1];
        elseif which_clust == 2
            this_col = [1 alpha_lin_time(max(1,n_counts)) alpha_lin_time(max(1,n_counts))];
        elseif which_clust == 3
            this_col = [alpha_lin_time(max(1,n_counts)) 1 ...
                alpha_lin_time(max(1,n_counts))];
        end
        %}
        
        if which_clust == 1
            this_col = [alpha_lin_time(max(1,n_counts)) ...
                alpha_lin_time(max(1,n_counts)) 1];
        elseif which_clust == 2
            this_col = [1 alpha_lin_time(max(1,n_counts)) alpha_lin_time(max(1,n_counts))];
        elseif which_clust == 3
            this_col = alpha_lin_time_yellow(:,max(1,n_counts))';
        end
        
        
        %alpha_temp = alpha_lin_time(max(1,n_counts));

        
        scatter3(locs(k,1),locs(k,2),locs(k,3),350,this_col,'filled');
        %alpha(p,alpha_temp);
        
        
    end
    xticklabels([])
    yticklabels([])
    zticklabels([])
    %xlabel('X')
    %ylabel('Y')
    %zlabel('Z')
    grid off
    
end


if 0
for tt = 1:nbins
    
    if sum(ch_counts(tt,:)) == 0
        continue
    end
    
    if i ~= 1
        additive_time = plot_times(group_idx(i,1))-plot_times(group_idx(1,1));
    else
        additive_time = 0;
    end
    
    
    %tt, i
    %round((bin_times(tt,1)-bin_times(1,1))/3600), round((bin_times(tt,1)-bin_times(1,1)+additive_time)/3600)
    
    
    fig = figure;
    set(gcf,'color','white');
    
    p = plotGIFTI(g);
    hold on
    
    % Get counts per channel in this time
    alpha_lin_time = linspace(1,0,max(ch_counts(tt,:)));
    alpha_lin_time_yellow = [linspace(1,0.9290,max(ch_counts(tt,:))); ...
        linspace(1,0.540,max(ch_counts(tt,:))); linspace(1,0.1250,max(ch_counts(tt,:)))];
    
    scatter3(locs(:,1),locs(:,2),locs(:,3),350,'k','linewidth',2);
    hold on
    %view(75.6 - tt*3,-4.2)
    if whichPt == 31
        view(75.6,-4.2)
    elseif whichPt == 9
        view(104.4,14.2)
    elseif whichPt == 17
        view(-196.3,1.2)
    elseif whichPt == 8
        view(-120,-11);
    end
    
    
    for k = 1:length(chs)
        
        % get the appropriate cluster
        for cl = 1:length(cl_locs)
            if ismember(k,cl_locs{cl}) == 1
                which_clust = cl;
            end
        end
        
        % Get appropriate color
        n_counts = ch_counts(tt,k);
        
        %{
        if which_clust == 1
            this_col = [alpha_lin_time(max(1,n_counts)) ...
                alpha_lin_time(max(1,n_counts)) 1];
        elseif which_clust == 2
            this_col = [1 alpha_lin_time(max(1,n_counts)) alpha_lin_time(max(1,n_counts))];
        elseif which_clust == 3
            this_col = [alpha_lin_time(max(1,n_counts)) 1 ...
                alpha_lin_time(max(1,n_counts))];
        end
        %}
        
        if which_clust == 1
            this_col = [alpha_lin_time(max(1,n_counts)) ...
                alpha_lin_time(max(1,n_counts)) 1];
        elseif which_clust == 2
            this_col = [1 alpha_lin_time(max(1,n_counts)) alpha_lin_time(max(1,n_counts))];
        elseif which_clust == 3
            this_col = alpha_lin_time_yellow(:,max(1,n_counts))';
        end
        
        
        %alpha_temp = alpha_lin_time(max(1,n_counts));

        
        scatter3(locs(k,1),locs(k,2),locs(k,3),350,this_col,'filled');
        %alpha(p,alpha_temp);
        
        
    end
    xticklabels([])
    yticklabels([])
    zticklabels([])
    %xlabel('X')
    %ylabel('Y')
    %zlabel('Z')
    grid off
    
    title(sprintf('Hour %d - %d',...
        round((bin_times(tt,1)-bin_times(1,1)+additive_time)/3600),...
        round((bin_times(tt,2)-bin_times(1,1)+additive_time)/3600)))
    set(gca,'fontsize',25)
    
    %set(gca,'visible','off')
    %title(sprintf('Hour %d',new_times(tt)),'fontsize',20);
    %{
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    %}
    
    F(tt) = getframe(fig);
    im = frame2im(F(tt));
    [imind,cm] = rgb2ind(im,256);
    
    if tt == 1
        imwrite(imind,cm,other_file_out,'gif', 'Loopcount',inf,'DelayTime',0.1);
    else
        imwrite(imind,cm,other_file_out,'gif','WriteMode','append','DelayTime',0.1);
    end

    close(fig)
    
    
end


end

end



end