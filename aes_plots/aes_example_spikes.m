function aes_example_spikes(pt,cluster)

whichPt = 31;

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
all_spikes(bad_idx) = [];
all_locs(bad_idx,:) = [];
idx(bad_idx) = [];
clusters = 1:k; clusters(bad_cluster) = [];
C(bad_cluster,:) = [];

colors = [0 0 1;1 0 0;0 1 0];
c_idx = zeros(size(idx,1),3);


end