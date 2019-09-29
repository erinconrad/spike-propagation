function predominant_cluster(pt,cluster,whichPts)

% Save file location
[~,~,scriptFolder,resultsFolder,~] = fileLocations;
destFolder = [resultsFolder,'clustering/stats/'];

if isempty(whichPts) == 1
    for i = 1:length(pt)
        if isempty(pt(i).seq_matrix) == 0
            if size(cluster(i).bad_cluster) < cluster(i).k
                whichPts = [whichPts,i];
            end
        end
    end
    
    % I should be doing all 20 of these patients
    if isequal(whichPts,[1,4,6,7,8,9,12,14,15,16,17,18,19,20,22,24,25,27,30,31]) == 0
     %   error('Warning, not doing correct patients!\n');
    end
end

%% Loop through patients
highest_rank_closest = [];
ranks_dists = cell(max(whichPts),1);
all_ranks_dists = [];
allAllDist = [];
allPredDist = [];
allNClust = [];

dist_most_freq = [];
dist_other = [];

for whichPt = whichPts
    
    fprintf('Doing %s\n',pt(whichPt).name);
    
    %% Get patient parameters
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
    
    %% Get cluster info
    
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
    
    %most popular cluster - note, I am identifying the element in the array
    %clusters that has the most popular cluster (instead of the cluster
    %number)
    popular = find(clusters == mode(idx)); 
    pop_ind = clusters == mode(idx);
    
    if size(C,1) <2, continue; end
    
    %{
    n_spikes_cluster = zeros(size(C,1),1);
    % get ranking of clusters by number of spikes
    for j = 1:size(C,1)
        n_spikes_cluster(j) = sum(idx == clusters(j));
    end
    [~,spike_ranking] = sort(n_spikes_cluster,'descend');
    %}
    
    
    % for each cluster, get the distance between its centroid and the
    % nearest SOZ
    soz_dist = zeros(size(C,1),1);
    for i = 1:size(C,1)
        soz_dist(i) = min(vecnorm(C(i,:) - locs(soz,:),2,2));    
    end
    
    [~,closest] = min(soz_dist);
    
    if popular == closest
        highest_rank_closest = [highest_rank_closest;1];
    else
        highest_rank_closest = [highest_rank_closest;0];
    end
    
    %{
    ranks_dists{whichPt} = [n_spikes_cluster,soz_dist];
    all_ranks_dists = [all_ranks_dists;n_spikes_cluster,soz_dist];
    
    
    % Distance from every electrode to its closest SOZ
    allLocs = zeros(size(locs,1),1);
    for i = 1:length(allLocs)
        allLocs(i) = min(vecnorm(locs(i,:) - locs(soz,:),2,2));
    end
    allAllDist =[allAllDist;mean(allLocs)];
    
    allPredDist = [allPredDist;soz_dist(popular)];
    %}
    allNClust = [allNClust;size(C,1)];
    
    dist_most_freq = [dist_most_freq;soz_dist(pop_ind)];
    dist_other = [dist_other;soz_dist(~pop_ind)];
 
end

% Number of patients for which most predominant is closest
num_closest = sum(highest_rank_closest);

%% Permutation test
if 0
% for each permutation, for each patient, randomly pick one of them
% clusters to identify as the most predominant cluster. Calculate the
% number of patiens where most popular cluster is closest. Do 1000 times
% and compare to actual result.

% Let's just say that the first cluster is "closest" and whenever I hit
% upon the first cluster, then it gets a point. I can do this because
% regardless of which one I define as closest, it's going to be the same
% result for a random permutation.

nb = 1e5; % I am doing 100,000 since p is close to 0.05
num_pop_closest = zeros(nb,1);
for ib = 1:nb
    
    % initialize number of patients where most popular cluster is closest
    num_pop_closest(ib) = 0;
    
    % Loop through all patients
    for i = 1:length(allNClust)
        
        % Get number of clusters for that patient
        n_clust = allNClust(i);
        
        % Pick random one
        x = randi(n_clust);
        
        % If x == 1, then I will say it is closest
        if x == 1
            num_pop_closest(ib) = num_pop_closest(ib) + 1;
        end
    end
end

% Find the number of times the permutation number is higher than the real
% number
num_pop_closest = sort(num_pop_closest);
diff_closest = num_pop_closest - num_closest;
allLarger = find(diff_closest > 0);
firstLarger = allLarger(1);

% The p value is then this number divided by the number of permutations
p = (nb - firstLarger+1)/(nb+1);

figure
histogram(num_pop_closest)
hold on
plot([num_closest num_closest],...
    get(gca,'ylim'));
text(sum(highest_rank_closest)+0.3,200,sprintf('%1.3f',p))
end


%% Second way - assume bernoulli random variable

% Expected value of number of patients where predominant cluster is closest
exp_num_closest = sum(1./allNClust);

% Variances add
var_num_closest = sum(1./(allNClust).*(1-1./allNClust));

%z = (num_closest-exp_num_closest)/sqrt(var_num_closest);

p = normcdf(num_closest,exp_num_closest,sqrt(var_num_closest),'upper')


%% Make plot
col1 = [0, 0.4470, 0.7410];
col2 = [0.8500 0.3250 0.0980];
figure
n = length(dist_most_freq);
xdata1 = ones(n,1).*rand(n,1)/2;
scatter(xdata1,dist_most_freq,150,col1,'filled')
hold on
plot([min(xdata1)-0.1 max(xdata1)+0.1],[median(dist_most_freq) median(dist_most_freq)],...
    'color',col1,'linewidth',3)
n2 = length(dist_other);
xdata2 = ones(n2,1).*rand(n2,1)/2+ones(n2,1);
scatter(xdata2,dist_other,150,col2,'filled')
plot([min(xdata2)-0.1 max(xdata2)+0.1],[median(dist_other) median(dist_other)],...
    'color',col2,'linewidth',3)
plot([(min(xdata1)+max(xdata1))/2,(min(xdata2)+max(xdata2))/2],...
    [max([dist_most_freq;dist_other])+10 max([dist_most_freq;dist_other])+10],...
    'k','linewidth',2)
text(((min(xdata1)+max(xdata1))/2+(min(xdata2)+max(xdata2))/2)/2-0.1,...
    max([dist_most_freq;dist_other])+15,'ns','fontsize',25)
ylim([-10 max([dist_most_freq;dist_other])+20]);
xlim([-0.2 1.7])
xticks([(min(xdata1)+max(xdata1))/2,(min(xdata2)+max(xdata2))/2]);
xticklabels({'Most frequent cluster','Other clusters'})
ylabel('Distance from seizure onset zone (mm)')
set(gca,'fontsize',20)

end