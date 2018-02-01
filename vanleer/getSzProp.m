function propSzInCluster = getSzProp(cluster_assignment,szOrNot,nclusters)

%{

This function calculates the proportion of spikes in each cluster that are
ictal

%}

propSzInCluster = zeros(10,1);

% loop through the clusters
for k = 1:nclusters
    
   % an array of length n_spikes with 1s where the spikes are in cluster k 
   temp = (cluster_assignment == k);
   
   % the total number in cluster k is the sum of the 1s
   n_incluster = sum(temp);
   
   % szOrNot(temp == 1) is an array with length equal to the number of spikes in
   % cluster k, with 1s at the spikes that are in seizures. By summing it
   % together, we get the number of spikes in that cluster that are in
   % seizures
   n_sz_incluster = sum(szOrNot(temp == 1));
   
   % we divide the number of spikes in that cluster that are in seizures by
   % the total number in the cluster to get the proportion of spikes in
   % that cluster that are in seizure
   propSzInCluster(k) = n_sz_incluster/n_incluster;
    
end


end