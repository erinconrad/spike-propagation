function clustersOverTime(pt,whichPts)

%% To do
%{
1) Think of removing sequences with mostly ties
2) Run over longer times
3) Run on peds data
4) Check out other clustering options


%}

%% Parameters

n_clusters = ones(20,1)*3;
%n_cluster = 


doPlots = 1;

% Save file location
[~,~,~,resultsFolder,~] = fileLocations;


for whichPt = whichPts

xyChan = pt(whichPt).electrodeData.locs;
saveFolder = [resultsFolder,'plots/',pt(whichPt).name,'/','clusters/'];
mkdir(saveFolder)

%% Get all sequences
[all_seq_cat,all_times] = divideIntoSzChunks(pt,whichPt);

all_seq_cat_old = all_seq_cat;
keep = ones(size(all_seq_cat,2),1);

%% Remove sequences with too many ties??
for s = 1:size(all_seq_cat_old,2)
   curr_seq = all_seq_cat_old(:,s);
   nonans = curr_seq(~isnan(curr_seq));
   norepeats = unique(nonans);
   if length(norepeats) < 0.5*length(nonans)
       keep(s) = 0;
   end
end

all_seq_cat(:,keep==0) = [];
all_times(:,keep == 0) = [];

fprintf('%s had %d sequences (%1.2f of all sequences) deleted for having >50 percent ties\n',...
    pt(whichPt).name,sum(keep == 0),sum(keep == 0)/length(keep));

%% Get all of the vectors
[all_vecs,early,late] = (getVectors2(all_seq_cat,pt(whichPt).electrodeData));
unit_vecs = all_vecs./vecnorm(all_vecs,2,2);

%% Get first channels
[~,firstChs] = min(all_seq_cat,[],1);
firstChs = xyChan(firstChs,2:4);

%% Cluster sequences by first channel and direction vector
final_vecs = unit_vecs;

for i = 1:30
[idx_all{i},C_all{i},sumd_all{i},D_all{i}] = ...
    kmeans([firstChs,final_vecs],n_clusters(whichPt));
metric(i) = sum(sumd_all{i});
end

% Take the results of the clustering algorithm that worked the best
[~,minidx] = min(metric);
idx = idx_all{minidx};
C = C_all{minidx};
D = D_all{minidx};

%% Get representative sequences
for i = 1:n_clusters(whichPt)
    [sortedD,I] = sort(D(:,i));
    rep_seq{i} = all_seq_cat(:,I(1:10));
    info(i).outputFile = [saveFolder,'cluster_',sprintf('%d',i),'.gif'];
    info(i).cluster = i;
    info(i).name = pt(whichPt).name;
end

%% Plot representative sequences
for i = 1:n_clusters
    movieSeqs(rep_seq{i},xyChan(:,2:4),info(i));
end

%% Assign each sequence a color based on what cluster index it it
colors = [1 0 0;0 1 0;0 0 1; 0.5 0.5 1; 1 0.5 0.5; 0.5 1 0.5];
c_idx = zeros(size(idx,1),3);
for i = 1:length(idx)
   c_idx(i,:) = colors(idx(i),:); 
end


%% Do plot
figure
set(gcf,'Position',[50 100 1200 1200])

% Plot of x, y, z coordinates of starting position over time
subplot(4,1,1)
toAdd = 0;
marker = {'x','o','>'};
for i = 1:3
scatter(all_times/3600,firstChs(:,1)+repmat(toAdd,size(firstChs,1),1),20,c_idx,marker{i})
hold on
if i ~=3
    toAdd = toAdd + quantile(firstChs(:,i),0.9) - quantile(firstChs(:,i+1),0.1);
end
end
for j = 1:length(pt(whichPt).sz)
   yl = ylim; 
   szOnset = pt(whichPt).sz(j).onset;
   sz = plot([szOnset szOnset]/3600,yl,'k','LineWidth',2);
end
set(gca,'ytick',[]);
title(sprintf('X, y, z coordinates of spike leader for %s',pt(whichPt).name));


% Plot of x, y, z coordinates of unit vector over time
subplot(4,1,2)
toAdd = 0;
marker = {'x','o','>'};
for i = 1:3
scatter(all_times/3600,final_vecs(:,1)+repmat(toAdd,size(final_vecs,1),1),20,c_idx,marker{i})
hold on
if i ~=3
    toAdd = toAdd + quantile(final_vecs(:,i),0.99) - quantile(final_vecs(:,i+1),0.01); 
end
end
set(gca,'ytick',[]);
for j = 1:length(pt(whichPt).sz)
   yl = ylim; 
   szOnset = pt(whichPt).sz(j).onset;
   sz = plot([szOnset szOnset]/3600,yl,'k','LineWidth',2);
end
title(sprintf('X, y, z coordinates of propagation vector for %s',pt(whichPt).name));

% Plot cluster identities
subplot(4,1,3)
scatter(all_times/3600,idx,20,c_idx);
hold on
for j = 1:length(pt(whichPt).sz)
   yl = ylim; 
   szOnset = pt(whichPt).sz(j).onset;
   sz = plot([szOnset szOnset]/3600,yl,'k','LineWidth',2);
end
title(sprintf('Cluster identities for %s',pt(whichPt).name));


% Plot locations of centroids 
subplot(4,1,4)
scatter3(xyChan(:,2),xyChan(:,3),xyChan(:,4),60,'k');
hold on

for k = 1:size(C,1)
    scatter3(C(k,1),C(k,2),C(k,3),60,colors(k,:),'filled');
    plot3([C(k,1) C(k,1) + C(k,4)],...
        [C(k,2) C(k,2) + C(k,5)],...
        [C(k,3) C(k,3) + C(k,6)],'k','LineWidth',2)
end
%}
title(sprintf('Spike leader and propagation vectors for %s',pt(whichPt).name));

saveas(gcf,[saveFolder,pt(whichPt).name,'cluster.png']);
close(gcf)

%}

end