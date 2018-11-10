function prettyClusterDist(pt,whichPt)

locs = pt(whichPt).electrodeData.locs(:,2:4);
cluster_vec = pt(whichPt).cluster.cluster_vec;
idx = pt(whichPt).cluster.idx;
k = pt(whichPt).cluster.k;

%% Assign each sequence a color based on what cluster index it it
colors = [0 0 1;1 0 0;0 1 0; 0.5 0.5 1; 1 0.5 0.5; 0.5 1 0.5; 0.4 0.7 0.4];
c_idx = zeros(size(idx,1),3);
for i = 1:length(idx)
   c_idx(i,:) = colors(idx(i),:); 
end


%% plots
figure
scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
hold on

scatter3(cluster_vec(:,1),cluster_vec(:,2),cluster_vec(:,3),...
    100,c_idx,'filled');


mean_vec = zeros(k,3);
for i = 1:k
   mean_vec(i,:) = mean([cluster_vec(idx==i,4),...
       cluster_vec(idx==i,5),cluster_vec(idx==i,6)],1);
end
    

xticklabels([])
yticklabels([])
zticklabels([])

view([0.7 0.2 0.2])


figure
scatter3(cluster_vec(:,4),cluster_vec(:,5),cluster_vec(:,6),...
    100,c_idx,'filled');

figure


end