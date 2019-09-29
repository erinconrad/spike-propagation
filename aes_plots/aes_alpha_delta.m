function aes_alpha_delta(pt,cluster,power,whichPt)


[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
p1 = genpath(scriptFolder);
addpath(p1);
destFolder = [resultsFolder,'alphaDelta/plots/'];
mkdir(destFolder);

    
szTimes = pt(whichPt).newSzTimes;
locs = pt(whichPt).electrodeData.locs(:,2:4);
soz = locs(pt(whichPt).newSOZChs,:);


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
all_spikes(bad_idx) = [];
all_locs(bad_idx,:) = [];
idx(bad_idx) = [];
clusters = 1:k; clusters(bad_cluster) = [];

% Get most popular cluster
popular = mode(idx);

bin_times = pt(whichPt).runTimes;
prop_pop = zeros(size(bin_times,1),1);
locs_bin = zeros(size(bin_times,1),3);
soz_dist_bin = zeros(size(bin_times,1),1);
num_spikes = zeros(size(bin_times,1),1);
SD_bin = zeros(size(bin_times,1),1);
std_z_bin = zeros(size(bin_times,1),1);
prop_pop_chunk = zeros(size(bin_times,1),1);
distNeeded = zeros(size(bin_times,1));
% Run through bin times and get proportion of spikes in most popular
% cluster for that bin.
for i = 1:size(bin_times,1)

    % Cluster identities of spikes in between those times
    whichClust = idx(all_times_all > bin_times(i,1) & ...
        all_times_all < bin_times(i,2));

    prop_pop(i) = sum(whichClust == popular)/length(whichClust);


    % Get spike locs
    whichSpikes = find(all_times_all > bin_times(i,1) & ...
        all_times_all < bin_times(i,2));
    locs_bin(i,:) = mean(all_locs(whichSpikes,:),1);
    num_spikes(i) = length(whichSpikes);



    % Get mean distance from spike to nearest SOZ
    soz_dist = zeros(length(whichSpikes),1);

    % Loop through all spikes in bin
    for j = 1:length(whichSpikes)

        % For each spike, get the distance between the spike and its
        % nearest SOZ
        soz_dist(j) = min(vecnorm(locs(all_spikes(whichSpikes(j)),:) - ...
            soz,2,2)); 
    end
    % Average that distance over all spikes in the bin
    soz_dist_bin(i) = mean(soz_dist);


    % Now get a measure of spatial dispersion for the bin
    n_spikes = length(whichSpikes); % number of spikes
    locs_sp = all_locs(whichSpikes,:); % location of every spike in the bin

    % Standard distance
    if size(locs_sp,1) < 20
        SD_bin(i) = nan;
    else
        SD_bin(i) = standardDistance(locs_sp);
    end
    std_z_bin(i) = std(locs_sp(:,3));

    % How many spikes are in most popular cluster for that chunk?
    prop_pop_chunk(i) = sum(whichClust == mode(whichClust))/length(whichClust);

    % distance from median to 60 percent of spikes in the chunk
    if size(locs_sp,1) < 10
        distNeeded(i) = nan;
    else
        distNeeded(i) = distToCaptureMost(locs_sp);
    end

end


%[badChNums,badChNamesOut] = getBadChs(pt,whichPt);
%mean_ad = mean(power(whichPt).ad_rat,1)';
all_ad = power(whichPt).alpha./power(whichPt).delta;
mean_ad = nanmean(all_ad,1);
%{
mean_ad = mean(power(whichPt).ad_rat(...
    ~ismember(1:length(pt(whichPt).channels),badChNums),:),1)';
%}



colors = [0 0 1;1 0 0;0 1 0];
c_idx = zeros(size(idx,1),3);



%% Prep for plot

plot_thing = all_locs;
plot_times = all_times_all;
sparse_time = plot_times-min(plot_times);

window = 3600;

for i = 1:length(clusters)
    clust{i} = plot_times(idx == clusters(i));
end

[sum_c,sum_times] = movingSumCounts(clust,plot_times,window);
totalSum = zeros(1,size(sum_times,2));
for i = 1:length(clusters)
    totalSum = totalSum + sum_c(i,:);
end
prop_c = sum_c./totalSum;

%% Plot

figure
set(gcf,'Position',[50 100 1200 500])
[ha, ~] = tight_subplot(2, 1, [.1 .01],[.08 .08],[.04 .01]);


axes(ha(1));
pl = zeros(length(clusters),1);
for i = 1:length(clusters)
    pl(i)= plot((sum_times-min(plot_times))/3600,prop_c(i,:),...
        'color',colors((i),:),'LineWidth',2);
hold on
end

possibleText = {'Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5','Cluster 6'};
textLeg = possibleText(1:size(C,1));

xlim([sparse_time(1)/3600-1 sparse_time(end)/3600+1])
set(gca,'xtick',[]);
set(gca,'FontSize',25);
title('Proportion of spikes in given cluster (moving average)')

axes(ha(2));
plot((power(whichPt).times-min(plot_times))/3600,nanmean(all_ad,1),...
    'k','LineWidth',2);
xlim([(plot_times(1)-min(plot_times))/3600-1 (plot_times(end)-min(plot_times))/3600+1])
hold on

% title(sprintf('Alpha-delta power ratio averaged across all electrodes'))

%% Re-label xticks
day0 = '01/01/00';
time0 = '10:58:42';
run_start_time = datetime([day0,' ',time0]);

% get seconds between xlim start and start_time
xlim_time_start = run_start_time + seconds(plot_times(1));

% get the datetime of the new rounded up xlim start time
start_round = xlim_time_start + hours(2) - ...
    minutes(minute(xlim_time_start)) - ...
    seconds(second(xlim_time_start));

% get the position (the time in hours)
xpos_start = seconds(start_round- xlim_time_start)/3600;

% now do the same thing for xlim end
xlim_time_end = run_start_time + seconds(plot_times(end));
end_round = xlim_time_end - hours(6)-...
    minutes(minute(xlim_time_end)) - ...
    seconds(second(xlim_time_end));
xpos_end = seconds(end_round- xlim_time_start)/3600;

% get a linear spacing of xtick positions
x_pos = linspace(xpos_start,xpos_end,13);

% get a linear spacing of times
xlim_times = linspace(start_round,end_round,13);

% turn these into strings
xlim_days = xlim_times.Day - min(xlim_times.Day) + 1;
xlim_hours = xlim_times.Hour;
xlim_str = {};
for i = 1:length(xlim_times)
    if xlim_hours(i) == 0
        new_hour = 12;
        pm = 'am';
    elseif xlim_hours(i) ==12
        new_hour = 12;
        pm = 'pm';
    elseif xlim_hours(i) >= 13
        new_hour = xlim_hours(i) - 12;
        pm = 'pm';
    else
        new_hour = xlim_hours(i);
        pm = 'am';
    end

    xlim_str{i} = sprintf('%d %s',new_hour,pm);

end
xticks(x_pos)
xticklabels(xlim_str)
title('Alpha/delta ratio')


set(gca,'FontSize',25);
legend([pl],[textLeg],'Position',[0.1 0.7 0.1 0.1])




end