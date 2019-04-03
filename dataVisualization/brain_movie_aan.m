function brain_movie_aan(pt,cluster,whichPt)


%% Parameters
nbins = 49;
delay = 0.2;

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;

outputFolder = [resultsFolder,'pretty_plots/'];
file_out = [outputFolder,'aan.gif'];
file_im_out = [outputFolder,'aan_im'];




%% Get data

 % Get patient parameters
locs = pt(whichPt).electrodeData.locs(:,2:4);
chs = 1:size(locs,1);
szTimes = pt(whichPt).newSzTimes;
soz = pt(whichPt).newSOZChs;


%% Load brain file
%{
brain_file = '/Users/erinconrad/Desktop/residency stuff/R25/actual work/data/brains/Study029/lh.pial.native.obj';
[V,F3,F4]=loadawobj(brain_file);

A = makeNewElecData(pt,whichPt);
%offset = [-0.4690 -50.8743 9.6682];
offset = [-0.4690 -50.8743 150]; 
locs = A*locs - offset;
%}

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


%% Bin times - number in each cluster over time
plot_times = all_times_all;
plot_thing = all_locs;
[Y,E] = discretize(plot_times,nbins);
new_times = round(E(1:end)/3600);
new_times = new_times-min(new_times) + 1;
new_counts = zeros(nbins,length(clusters));
for bb = 1:nbins
    for k = 1:length(clusters)
        new_counts(bb,k) = sum(Y==bb & idx==clusters(k));
    end
end
new_prop = new_counts./sum(new_counts,2);

%% Bin times - number in each channel over time
plot_times = all_times_all;
[Y,E] = discretize(plot_times,nbins);
ch_counts = zeros(nbins,length(chs));
for bb = 1:nbins
    for k = 1:length(chs)
        ch_counts(bb,k) = sum(Y==bb & all_spikes==k);
    end
end

%cols = [0 0 1;1 0 0;0 1 0];

alpha_lin = linspace(0.1,1,max(max(ch_counts)));


%% plot figure showing locations of spikes and clusters
fig = figure;

for k = 1:length(chs)
    % get the appropriate cluster
    for cl = 1:length(cl_locs)
        if ismember(k,cl_locs{cl}) == 1
            which_clust = cl;
        end
    end

        if which_clust == 1
            this_col = [0 0 1];
        elseif which_clust == 2
            this_col = [1 0 0];
        elseif which_clust == 3
            this_col = [0 1 0];
        end

        scatter3(locs(k,1),locs(k,2),locs(k,3),200,this_col,'filled');
        hold on

end
scatter3(locs(:,1),locs(:,2),locs(:,3),200,'k','linewidth',2);
view(75.6,-4.2)
xticklabels([])
yticklabels([])
zticklabels([])
grid off
set(gca,'visible','off')
pause
print(fig,file_im_out,'-depsc');
close(fig)


%% Plot movie


for tt = 1:nbins
    
    fig = figure;
   % set(fig,'Position',[100 100 
    set(gcf,'color','white');
    
    %{
    % Plot the brain
    patch('Vertices',V','Faces',F3','FaceColor',[0.9 0.75 0.75],'EdgeColor','None');
    hold on
    camlight(gca,-80,-10);
    lighting(gca,'gouraud');
    %}
    
    % Get counts per channel in this time
    alpha_lin_time = linspace(1,0,max(ch_counts(tt,:)));
    
    scatter3(locs(:,1),locs(:,2),locs(:,3),350,'k','linewidth',2);
    hold on
    view(75.6,-4.2)
    
    for k = 1:length(chs)
        
        % get the appropriate cluster
        for cl = 1:length(cl_locs)
            if ismember(k,cl_locs{cl}) == 1
                which_clust = cl;
            end
        end
        
        % Get appropriate color
        n_counts = ch_counts(tt,k);
        
        if which_clust == 1
            this_col = [alpha_lin_time(max(1,n_counts)) ...
                alpha_lin_time(max(1,n_counts)) 1];
        elseif which_clust == 2
            this_col = [1 alpha_lin_time(max(1,n_counts)) alpha_lin_time(max(1,n_counts))];
        elseif which_clust == 3
            this_col = [alpha_lin_time(max(1,n_counts)) 1 ...
                alpha_lin_time(max(1,n_counts))];
        end
        
        
        %alpha_temp = alpha_lin_time(max(1,n_counts));

        
        scatter3(locs(k,1),locs(k,2),locs(k,3),350,this_col,'filled');
        %alpha(p,alpha_temp);
        
        
    end
    xticklabels([])
    yticklabels([])
    zticklabels([])
    grid off
    set(gca,'visible','off')
    title(sprintf('Hour %d',new_times(tt)),'fontsize',20);
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
    
    F(tt) = getframe(fig);
    im = frame2im(F(tt));
    [imind,cm] = rgb2ind(im,256);
    
    if tt == 1
        imwrite(imind,cm,file_out,'gif', 'Loopcount',inf,'DelayTime',delay);
    else
        imwrite(imind,cm,file_out,'gif','WriteMode','append','DelayTime',delay);
    end

    close(fig)
    
end


end