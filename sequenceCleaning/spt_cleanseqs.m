function [clueGraph,seq_track] = spt_cleanseqs(clueGraph)

% 'Noisy' sequences are detected using Kmeans clustering 
% and the weighted degree distribution. 

%Perform kmeans clustering
seq_track  = 1:size(clueGraph,2);                 %track removed sequences
in_deg     = sum(clueGraph,1);     
out_deg    = sum(clueGraph,2);     
degree     = in_deg + out_deg';    
K_data     = [in_deg',out_deg];                   %kmeans inputs
[idx,C]    = kmeans(K_data,3);                    %perform kmeans

%Which idx contains the noisy seqs?
grp1     = mean(degree(find(idx==1)));
grp2     = mean(degree(find(idx==2)));
grp3     = mean(degree(find(idx==3)));
grps     = [grp1,grp2,grp3];
[~,toss] = min(grps);

%Toss noisy sequences from clueGraph
toss_idx                     = find(idx==toss)';
clueGraph(:,toss_idx)        = [];
clueGraph(toss_idx,:)        = [];
seq_track(toss_idx)          = [];
