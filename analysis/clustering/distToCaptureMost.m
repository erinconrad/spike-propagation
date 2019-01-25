function distNeeded = distToCaptureMost(locs)

thresh = 0.6;

n = size(locs,1);
dist_sp = zeros(n,1);
n_thresh = round(thresh*n);

%% Get location of the median
median_loc = median(locs,1);


%% Loop through spikes and get distance of each spike from median
for i = 1:size(locs,1)
    dist_sp(i) = norm(locs(i,:) - median_loc);    
end

%% Sort the distances
sorted_distances = sort(dist_sp);


%% Distance needed is distance of 80th percentile
if n_thresh == 0
    distNeeded = nan;
else
    distNeeded = sorted_distances(n_thresh);
end

end