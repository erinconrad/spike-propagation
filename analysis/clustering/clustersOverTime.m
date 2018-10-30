function clustersOverTime(pt,whichPts)

%% To do
%{
1) Think of other ways to validate
2) Run over longer times
3) Run on peds data
4) Check out other clustering options
5) Come up with some statistical tests
      - does the sequence cluster vary over time (over 30 minute chunks,
      perhaps) (the answer will be yes, kind of exciting)
      - are seizures more likely during certain cluster distributions?
      - does the cluster distribution change prior to seizures?
6) run on ictal data
7) Come up with nicer visualization for poster and to show Eric


%}

%% FYIs
%{

- I clip ten minutes before to ten minutes after the seizure in divideIntoSzChunks

%}

%% Parameters
window = 3600;
n_clusters = ones(20,1)*3;


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

%% Cluster sequences by first channel and direction vecto
final_vecs = unit_vecs;

cluster_approach = 1;
if cluster_approach == 1

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
else
    [clustCent,data2cluster,cluster2dataCell] = MeanShiftCluster([firstChs,final_vecs]',50);
    C = clustCent;
    idx = data2cluster;
    n_clusters(whichPt) = length(unique(idx));
end

%% Get representative sequences

for i = 1:n_clusters(whichPt)
    [sortedD,I] = sort(D(:,i));
    rep_seq{i} = all_seq_cat(:,I(1:10));
    info(i).outputFile = [saveFolder,'cluster_',sprintf('%d',i),'.gif'];
    info(i).cluster = i;
    info(i).name = pt(whichPt).name;
end

%% Plot representative sequences
%{
for i = 1:n_clusters
    movieSeqs(rep_seq{i},xyChan(:,2:4),info(i));
end
%}
%}

%% Assign each sequence a color based on what cluster index it it
colors = [0 0 1;1 0 0;0 1 0; 0.5 0.5 1; 1 0.5 0.5; 0.5 1 0.5; 0.4 0.7 0.4];
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
    toAdd = toAdd + quantile(final_vecs(:,i),0.9999) - quantile(final_vecs(:,i+1),0.0001); 
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

subplot(5,1,3)
scatter(all_times/3600,idx,10,c_idx,'filled');
hold on
for j = 1:length(pt(whichPt).sz)
   yl = ylim; 
   szOnset = pt(whichPt).sz(j).onset;
   sz = plot([szOnset szOnset]/3600,yl,'k','LineWidth',2);
end
title(sprintf('Cluster identities for %s',pt(whichPt).name));
%}

%% Plot proportion of sequences in a given cluster over a moving window
% Moving sum
for i = 1:n_clusters(whichPt)
clust{i} = all_times(idx == i);
end

[sum_c,sum_times] = movingSumCounts(clust,all_times,window);

totalSum = zeros(1,size(sum_times,2));
for i = 1:n_clusters(whichPt)
    totalSum = totalSum + sum_c(i,:);
end

subplot(4,1,3)


for i = 1:n_clusters(whichPt)
   pl(i)= plot(sum_times/3600,sum_c(i,:)./totalSum,'color',colors(i,:),'LineWidth',2);
%plot(chunk_times,n_clusters_chunk(:,i),'color',colors(i,:),'LineWidth',2);
hold on
end

for j = 1:length(pt(whichPt).sz)
   yl = ylim; 
   szOnset = pt(whichPt).sz(j).onset;
   sz = plot([szOnset szOnset]/3600,yl,'k','LineWidth',2);
end

%legend([pl(1) pl(2) pl(3)],{'Cluster 1','Cluster 2','Cluster 3'});
title(sprintf(['Proportion of sequences in given cluster, moving'...
    ' average %d s, %s'],...
    window,pt(whichPt).name));

%}


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

title(sprintf('Spike leader and propagation vectors for %s',pt(whichPt).name));

saveas(gcf,[saveFolder,pt(whichPt).name,'cluster.png']);
%close(gcf)

%}

%}

%% Statistical tests

% Get all seizure times
szAll = zeros(length(pt(whichPt).sz),1);
for j = 1:length(pt(whichPt).sz)
   szAll(j) = pt(whichPt).sz(j).onset;
    
end

%% #1 does cluster distribution change over 60 minute chunks?

test_t = 3600;

% Determine most popular cluster
pop_c = mode(idx);

% Divide run into 60 minute chunks
n_chunks = ceil((all_times(end) - all_times(1))/test_t);

prop_pop = zeros(n_chunks,1);
times_pop = zeros(n_chunks,1);
which_chunk = zeros(length(all_times),1);
num_cluster = zeros(n_chunks,n_clusters(whichPt));
sz_chunk = zeros(n_chunks,1);
most_num = zeros(n_chunks,1);

for i = 1:n_chunks
   curr_times = [all_times(1)+(i-1)*test_t, min(all_times(1) + i*test_t,all_times(end))];
   curr_seqs = find(all_times >= curr_times(1) & all_times <= curr_times(2));
   prop_pop(i) = sum(idx(curr_seqs) == pop_c)/length(curr_seqs);
   which_chunk(curr_seqs) = i;
   times_pop(i) = curr_times(2);
   for j = 1:size(num_cluster,2)
      num_cluster(i,j) = sum(idx(curr_seqs) == j); 
       
   end
   
   if any(szAll >= curr_times(1) & szAll <= curr_times(2))
       sz_chunk(i) = 1;
   end
   
   most_num(i) = mode(idx(curr_seqs));
end



% Do an chi-squared to test if the 
% cluster changes across the 60 minute chunks
[tbl_1,chi2_1,p_1,labels_1] = crosstab(which_chunk,idx);


%% #2 Are 60 minute chunks containing seizures more likely to have certain cluster distributions
[tbl_2,chi2_2,p_2,labels_2] = crosstab(sz_chunk,most_num);

end

