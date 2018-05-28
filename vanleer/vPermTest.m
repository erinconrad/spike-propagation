function stats =  vPermTest(cluster_assignment,szOrNot,nclusters,nboot)

%{

This function takes the cluster assignments and ictal/preictal assignments
for all spikes and does a permutation test to answer whether there is a
significantly different proportion in each cluster of spikes that are ictal
versus preictal

%}


alpha = 0.01;

% This is how many of the smallest and largest proportions from the
% permutations to take as the significance confidence interval. Basically,
% anything more extreme than these will be considered significant
nOutliers = nboot*alpha/2;

propSzInClusterBoot = zeros(nboot,nclusters);

% Do it once for real
% Get proportion of spikes in each cluster that occur during a seizure
propSzInClusterReal = (getSzProp(cluster_assignment,szOrNot,nclusters))';

% Loop through the number of permutations
for i = 1:nboot
    if mod(i,1e3) == 0
        fprintf('On permutation %d\n',i);
    end
    
    % randomize the identities of seizure/not seizure
    p = randperm(length(szOrNot));
    tempSzAssignments = szOrNot(p);
    
    % get the proportion of spikes in each cluster that are ictal for this
    % fake data
    propSzInCluster = getSzProp(cluster_assignment,tempSzAssignments,nclusters);
    propSzInClusterBoot(i,:) = propSzInCluster;
end

%% Get some stats on the bootstrap

% the mean of the bootstrap data (I would expect this to be about equal
% proportions across all clusters)
meanProp = mean(propSzInClusterBoot,1); 

% Sort by proportion
propSort = sort(propSzInClusterBoot,1);


limsProp = [propSort(nOutliers,:);propSort(end-nOutliers,:)];

stats.real = propSzInClusterReal;
stats.meanBoot = meanProp;
stats.limsProp = limsProp;
end