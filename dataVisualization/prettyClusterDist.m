function prettyClusterDist(pt,whichPt)

%% Get paths and load seizure info and channel info
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
p1 = genpath(scriptFolder);
addpath(p1);

destFolder = [resultsFolder,'pretty_plots/Fig2/'];

s=20;
sizey = 400;

locs = pt(whichPt).electrodeData.locs(:,2:4);
xyChan = pt(whichPt).electrodeData.locs;
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
set(gcf,'position',[200 200 1000 800]);

scatter3(locs(:,1),locs(:,2),locs(:,3),sizey,'k','linewidth',3);
hold on
scatter3(cluster_vec(:,1),cluster_vec(:,2),cluster_vec(:,3),...
    sizey,c_idx,'filled');

hold on

%% Re-plot a couple for the purpose of making the legend
t = find(idx == 1);
t = t(1);
pl1 = scatter3(cluster_vec(t,1),cluster_vec(t,2),cluster_vec(t,3),...
    sizey,c_idx(t,:),'filled');

t = find(idx == 2);
t = t(1);
pl2 = scatter3(cluster_vec(t,1),cluster_vec(t,2),cluster_vec(t,3),...
    sizey,c_idx(t,:),'filled');

t = find(idx == 3);
t = t(1);
pl3 = scatter3(cluster_vec(t,1),cluster_vec(t,2),cluster_vec(t,3),...
    sizey,c_idx(t,:),'filled');





mean_vec = zeros(k,3);
for i = 1:k
   mean_vec(i,:) = mean([cluster_vec(idx==i,4),...
       cluster_vec(idx==i,5),cluster_vec(idx==i,6)],1);
end
    
%% Add vectors
%{
didCh = zeros(size(locs,1),1);
for i = 1:length(idx)
    chLoc = cluster_vec(i,1:3);
    for j = 1:size(xyChan,1)
       if isequal(xyChan(j,2:4),chLoc) == 1
           if didCh(xyChan(j,1)) == 1
               continue
           end
           didCh(xyChan(j,1)) = 1;
           whichClust = idx(i);
           quiver3(xyChan(j,2),xyChan(j,3),xyChan(j,4),...
               s*mean_vec(whichClust,1),s*mean_vec(whichClust,2),s*mean_vec(whichClust,3),...
               'color',colors(whichClust,:),'linewidth',5,'maxheadsize',20);
           
       end
        
    end
    
end
%}

%% Plot SOZ channel
soz = pt(whichPt).newSOZChs;
scatter3(locs(soz,1),locs(soz,2),locs(soz,3),...
    100,'k','filled');

%{
legend([pl1,pl2,pl3],{'Cluster 1','Cluster 2','Cluster 3'},'FontSize',50,...
    'Location','southeast');
    %}

xticklabels([])
yticklabels([])
zticklabels([])

view([0.7 0.2 0.2])
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 1000/800*20 20];
print(gcf,[destFolder,'clustEx111'],'-depsc');
print(gcf,[destFolder,'clustEx111'],'-dpng');
eps2pdf([destFolder,'clustEx111.eps'])


end