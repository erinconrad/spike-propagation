function clustersOverTime(pt,whichPts)

%% Parameters
allSpikes = 1; % Instead of just lead spike, look at all spikes
clustOpt = 0; % plot elbow plot
doPlots = 1; % do main plots
doLongPlots = 1; % do long plots
removeTies = 1; %remove sequences containing too many ties
leadOnly = 1; %don't change
window = 3600;


%% Optimal cluster numbers
%n_clusters = ones(30,1)*3;
if allSpikes == 1
    n_clusters(3) = 4; %HUP68, 2 by silhouette
    n_clusters(4) = 2; %HUP70
    n_clusters(8) = 3; %HUP78
    n_clusters(9) = 3; %HUP080
    n_clusters(12) = 3; %HUP86
    n_clusters(17) = 2; %HUP106
    n_clusters(18) = 4; %HUP107
    n_clusters(19) = 3; %HUP111A
    n_clusters(20) = 3; %HUP116
    n_clusters(22) = 4; %Study16
    n_clusters(24) = 3; %Study19
    n_clusters(25) = 4; %Study20
    n_clusters(27) = 3; %Study22
    n_clusters(30) = 4; %Study28
    n_clusters(31) = 3; %Study29
else
    n_clusters(3) = 4; %HUP68
    n_clusters(4) = 5; %HUP70
    n_clusters(8) = 3; %HUP78
    n_clusters(9) = 4; %HUP080
    n_clusters(12) = 5; %HUP86
    n_clusters(17) = 3; %HUP106
    n_clusters(18) = 4; %HUP107
    n_clusters(19) = 3; %HUP111A
    n_clusters(20) = 5; %HUP116
    n_clusters(22) = 3; %Study16
    n_clusters(24) = 3; %Study19
    n_clusters(25) = 4; %Study20
    n_clusters(27) = 4; %Study22
    n_clusters(30) = 4; %Study28
end


% Save file location
[~,~,~,resultsFolder,~] = fileLocations;


for whichPt = whichPts
   
fprintf('Doing %s\n',pt(whichPt).name);
if isempty(pt(whichPt).electrodeData) == 1
    continue
end

xyChan = pt(whichPt).electrodeData.locs;
saveFolder = [resultsFolder,'cluster_validation_all/',pt(whichPt).name,'/'];
mkdir(saveFolder)

%% Get all sequences

[all_seq_cat,all_times,~,~] = divideIntoSzChunksGen(pt,whichPt);

all_seq_cat_old = all_seq_cat;

% Test that I removed all ictal and periictal sequences
%{
szTimes = zeros(length(pt(whichPt).sz),2);
for j = 1:length(pt(whichPt).sz)
    szTimes(j,:) = [pt(whichPt).sz(j).onset - 60 pt(whichPt).sz(j).offset + 60];
end

firstSp = min(all_seq_cat,[],1);
t = find(any(firstSp >= szTimes(:,1) & firstSp <= szTimes(:,2)));
%}


%% Remove sequences with too many ties??
keep = ones(size(all_seq_cat,2),1);
if removeTies == 1
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

fprintf(['%s had %d sequences (%1.2f of all sequences) deleted'...
    'for having >50 percent ties\n%d sequences remain\n'],...
    pt(whichPt).name,sum(keep == 0),sum(keep == 0)/length(keep),sum(keep==1));
end

%% Add these sequences to the pt Struct
pt(whichPt).cluster.all_seq_cat = all_seq_cat;
pt(whichPt).cluster.all_times = all_times;

%% Get all of the vectors
[all_vecs,early,late] = (getVectors2(all_seq_cat,pt(whichPt).electrodeData));
unit_vecs = all_vecs./vecnorm(all_vecs,2,2);

%% Remove instances of zero vectors
I = find(vecnorm(all_vecs,2,2)==0);
if isempty(I) == 0
    fprintf('warning, some zero length vectors\n');
    all_seq_cat(:,I) = [];
    all_times(:,I) = [];
    unit_vecs(I,:) = [];
    all_vecs(I,:) = [];
end

%% Get first channels
[~,firstChs] = min(all_seq_cat,[],1);
firstChs = xyChan(firstChs,2:4);
[~,lastChs] = max(all_seq_cat,[],1);
lastChs = xyChan(lastChs,2:4);
final_vecs = unit_vecs;
cluster_vec = [firstChs,final_vecs];

pt(whichPt).cluster.cluster_vec = cluster_vec;
pt(whichPt).cluster.firstChs = firstChs;
pt(whichPt).cluster.unit_vecs = unit_vecs;
pt(whichPt).cluster.k = n_clusters(whichPt);

%% Get all spikes
all_spikes = [];
all_times_all = [];
seq_index = [];
for i = 1:size(all_seq_cat,2)
    nonan = find(~isnan(all_seq_cat(:,i)));
    all_spikes = [all_spikes;nonan];
    all_times_all = [all_times_all;all_seq_cat(nonan,i)];
    seq_index = [seq_index;i*ones(length(nonan),1)];
end
all_locs = xyChan(all_spikes,2:4);
pt(whichPt).cluster.all_spikes = all_spikes;
pt(whichPt).cluster.all_locs = all_locs;
pt(whichPt).cluster.all_times_all = all_times_all;

%% Determine optimal number of clusters


% Elbow approach
if clustOpt == 1
SSE = zeros(10,1);
for k = 1:10
    
    SSE_temp = zeros(30,1);
    for j = 1:30
        if allSpikes == 1
            [idx_test,C_test] = ...
                kmeans(all_locs,k);
        elseif leadOnly == 0
            [idx_test,C_test] = ...
                kmeans(cluster_vec,k);
        else
            [idx_test,C_test] = ...
                kmeans(firstChs,k);
        end

        % Get SSE
        for i = 1:k
            if allSpikes == 1
                SSE_temp(j) = SSE_temp(j) + sum(sum((all_locs(idx_test == i,:) - ...
                   repmat(C_test(i,:),size(all_locs(idx_test == i,:),1),1)).^2));
            elseif leadOnly == 0
               SSE_temp(j) = SSE_temp(j) + sum(sum((cluster_vec(idx_test == i,:) - ...
                   repmat(C_test(i,:),size(cluster_vec(idx_test == i,:),1),1)).^2));
            else
                SSE_temp(j) = SSE_temp(j) + sum(sum((firstChs(idx_test == i,:) - ...
                   repmat(C_test(i,:),size(firstChs(idx_test == i,:),1),1)).^2));
            end
            
        end

    end
    SSE(k) = min(SSE_temp);
    
    %{
    % For Erin's education
    idx = sub2ind(size(C_test), repmat(idx_test, [1 p]), ...
        repmat(1:p, [length(idx_test) 1]));
    C = C_test(idx);
    SSE(k) = sum(sum((cluster_vec - C).^2));
    %}
    
end

figure
plot(1:10,SSE)
end

% Silhouette method
if 1 == 0
E = evalclusters(all_locs,'kmeans','silhouette','klist',[1:10]);
pt(whichPt).cluster.optimalKSilhouette = E.OptimalK; 
%gscatter(cluster_vec(:,1),cluster_vec(:,2),E.OptimalY,'rbg','xod')
end

%{
myfunc = @(X,K)(kmeans(X, K, 'emptyaction','singleton',...
    'replicate',5));
eva = evalclusters([firstChs,final_vecs],myfunc,'CalinskiHarabasz',...
    'klist',[1:6]);
plot(eva)
eva
%}

%E = evalclusters([firstChs,final_vecs],'kmeans','DaviesBouldin','klist',[1:6]);

% Gap method
if 1 == 0
eva = evalclusters([firstChs],'kmeans','gap','KList',[1:20]);
end

%% Do clustering algorithm

cluster_approach = 1;
if cluster_approach == 1
    % K means
    for i = 1:30
        if allSpikes == 1
        [idx_all{i},C_all{i},sumd_all{i},D_all{i}] = ...
        kmeans([all_locs],n_clusters(whichPt));
        elseif leadOnly == 0
    [idx_all{i},C_all{i},sumd_all{i},D_all{i}] = ...
        kmeans([firstChs,final_vecs],n_clusters(whichPt));
        else
     [idx_all{i},C_all{i},sumd_all{i},D_all{i}] = ...
        kmeans([firstChs],n_clusters(whichPt));      
        end
        
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

pt(whichPt).cluster.idx = idx;
pt(whichPt).cluster.C = C;
pt(whichPt).cluster.D = D;

%% Get representative sequences

for i = 1:n_clusters(whichPt)
    [sortedD,I] = sort(D(:,i));
    
    if allSpikes == 1
        % I want to find the sequence that the spike came from.
        whichSeqs = seq_index(I(1:12));
        rep_seq{i} = all_seq_cat(:,whichSeqs);
    else
        rep_seq{i} = all_seq_cat(:,I(1:12));
    end

    info(i).outputFile = [saveFolder,'cluster_',sprintf('%d',i),'.gif'];
    info(i).cluster = i;
    info(i).name = pt(whichPt).name;
end

pt(whichPt).cluster.rep_seq = rep_seq;

%% Plot representative sequences
if doLongPlots == 1
    for i = 1:n_clusters(whichPt)
        outputFile = [saveFolder,sprintf('seqs_cluster_%d',i),'.png'];
        showSpecificSequences(pt,whichPt,rep_seq{i},1,outputFile)
    end



    for i = 1:n_clusters(whichPt)
        movieSeqs(rep_seq{i},xyChan(:,2:4),info(i));
    end
end



%% Assign each sequence a color based on what cluster index it it
colors = [0 0 1;1 0 0;0 1 0; 0.5 0.5 1; 1 0.5 0.5; 0.5 1 0.5; 0.4 0.7 0.4];
c_idx = zeros(size(idx,1),3);
for i = 1:length(idx)
   c_idx(i,:) = colors(idx(i),:); 
end

if doPlots == 1
%% Do plot
figure
set(gcf,'Position',[50 100 1200 1200])

% Plot of x, y, z coordinates of starting position over time
subplot(3,1,1)
toAdd = 0;
%marker = {'x','o','>'};
ttext = {'x','y','z'};

if allSpikes == 1
    plot_thing = all_locs;
    plot_times = all_times_all;
else
    plot_thing = firstChs;
    plot_times = all_times;
end

for i = 1:3
scatter(plot_times/3600,plot_thing(:,i)+repmat(toAdd,size(plot_thing,1),1),20,c_idx)
hold on
text(plot_times(1)/3600-0.3,toAdd+median(plot_thing(:,i)),sprintf('%s',ttext{i}),'FontSize',30);
if i ~=3
    toAdd = toAdd + 10+(max(plot_thing(:,i)) - min(plot_thing(:,i+1)));%quantile(firstChs(:,i),0.95) - quantile(firstChs(:,i+1),0.05);
end
end
for j = 1:length(pt(whichPt).sz)
   yl = ylim; 
   szOnset = pt(whichPt).sz(j).onset;
   sz = plot([szOnset szOnset]/3600,yl,'k','LineWidth',2);
end
set(gca,'ytick',[]);
xlim([plot_times(1)/3600-1 plot_times(end)/3600+1])
title(sprintf('X, y, z coordinates of spike leader for %s',pt(whichPt).name));
set(gca,'FontSize',15);

if 1 == 0
% Plot of x, y, z coordinates of unit vector over time
subplot(4,1,2)
toAdd = 0;
%marker = {'x','o','>'};
for i = 1:3
    scatter(all_times/3600,final_vecs(:,i)+repmat(toAdd,size(final_vecs,1),1),20,c_idx)
    hold on
    
    text(all_times(1)/3600-0.3,toAdd,sprintf('%s',ttext{i}),'FontSize',30);
    if i ~=3
        toAdd = toAdd + 3;%quantile(final_vecs(:,i),0.9999) - quantile(final_vecs(:,i+1),0.0001); 
    end
end
set(gca,'ytick',[]);
for j = 1:length(pt(whichPt).sz)
   yl = ylim; 
   szOnset = pt(whichPt).sz(j).onset;
   sz = plot([szOnset szOnset]/3600,yl,'k','LineWidth',2);
end
xlim([all_times(1)/3600-1 all_times(end)/3600+1])
title(sprintf('X, y, z coordinates of propagation vector for %s',pt(whichPt).name));
set(gca,'FontSize',15);
end

% Plot cluster identities
%{
subplot(4,1,3)
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
clust{i} = plot_times(idx == i);
end

[sum_c,sum_times] = movingSumCounts(clust,plot_times,window);

totalSum = zeros(1,size(sum_times,2));
for i = 1:n_clusters(whichPt)
    totalSum = totalSum + sum_c(i,:);
end

subplot(3,1,2)

pl = zeros(n_clusters(whichPt),1);
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
xlim([plot_times(1)/3600-1 plot_times(end)/3600+1])

leg_text = {};
for i = 1:length(pl)
   leg_text = [leg_text,sprintf('Cluster %d',i)];
end

legend([pl],...
    leg_text,'Position',...
    [0.87 0.9 0.1 0.05]);
title(sprintf(['Proportion of sequences in given cluster, moving'...
    ' average %d s, %s'],...
    window,pt(whichPt).name));
set(gca,'FontSize',15);

%}


% Plot locations of centroids 


subplot(3,1,3)
scatter3(xyChan(:,2),xyChan(:,3),xyChan(:,4),60,'k');
hold on

for k = 1:size(C,1)
    scatter3(C(k,1),C(k,2),C(k,3),60,colors(k,:),'filled');
    if 1 == 0
    plot3([C(k,1) C(k,1) + 10*C(k,4)],...
        [C(k,2) C(k,2) + 10*C(k,5)],...
        [C(k,3) C(k,3) + 10*C(k,6)],'k','LineWidth',2)
    end
end

title(sprintf('Spike leader and propagation vectors for %s',pt(whichPt).name));
set(gca,'FontSize',15);
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'zticklabel',[])
%saveas(gcf,[saveFolder,pt(whichPt).name,'cluster.png']);
%close(gcf)

end
%}



%{

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

fprintf(['For %s, regarding whether 60 minute chunks\n have different cluster'...
    ' distributions,\n the p-value is %1.1e\n\n\n'],pt(whichPt).name,p_1);


%% Do test to compare cluster distributions between pre-ictal and inter-ictal data

% Get all the pre-ictal sequences
preIcRange = [-40*60,-10*60];
preIcClustIdx = [];
preIcIdx = [];
preIcTimes = zeros(length(pt(whichPt).sz),2);
for j = 1:length(pt(whichPt).sz)
    szTime = pt(whichPt).sz(j).onset; 
    preIcTime = szTime + preIcRange;
    preIcTimes(j,:) = preIcTime;
    preIcClustIdx = [preIcClustIdx;idx(all_times >= preIcTime(1) & ...
        all_times <= preIcTime(2))];
    preIcIdx = [preIcIdx;find(all_times >= preIcTime(1) & ...
        all_times <= preIcTime(2))'];
    
end

% The harder part: get random 30 minute chunks in the interictal period
% equivalent to the number of pre-ictal chunks
szTimes = zeros(length(pt(whichPt).sz),1);
for j= 1:length(pt(whichPt).sz)
    szTimes(j) = pt(whichPt).sz(j).onset; 
end

%{
n_chunks = 3*length(pt(whichPt).sz);
interIcTimes = zeros(n_chunks,2);
i_chunk = 1;
while i_chunk <= n_chunks
    t_1 = randi([round(all_times(1)),round(all_times(end))]);
    if any(abs(interIcTimes(:,1)-t_1)<=1*3600)
        continue;
    end
    
    if any(abs(szTimes-t_1) <= 3*3600)
        continue
    end
    
    if any(abs(t_1-szTimes) <= 1*3600)
        continue
    end
    
    interIcTimes(i_chunk,:) = [t_1 t_1+60*60];
    i_chunk = i_chunk + 1;
    
end
%}

% Alternate idea without randomness:
n_chunks = length(pt(whichPt).sz);
interIcTimes = [];
for i = 1:n_chunks
    
    % potential late_time is 4 hours before the pre-ictal period
   late_time = preIcTimes(i,1) - 3600*2;
   
   % potential early time is either the start of the run or an hour after
   % the last seizure
   if i ==1
       early_time = all_times(1);
   else
       early_time = szTimes(i-1) + 3600*1;
   end
   
   if early_time > late_time
       continue
   end
   
   interIcTimes = [interIcTimes;early_time late_time];

end
 % add another time
 late_time = all_times(end);
 early_time = szTimes(end) + 3600*1;
 if late_time>early_time
     interIcTimes = [interIcTimes;early_time late_time];
 end

interIcClustIdx = [];
interIcIdx = [];
for i = 1:size(interIcTimes,1)
    interIcClustIdx = [interIcClustIdx;idx(all_times >= interIcTimes(i,1) & ...
        all_times <= interIcTimes(i,2))];
    interIcIdx = [interIcIdx;find(all_times >= interIcTimes(i,1) & ...
        all_times <= interIcTimes(i,2))'];

    
end

[tbl_2,chi2_2,p_2,labels_2] = crosstab([ones(size(preIcClustIdx));...
    2*ones(size(interIcClustIdx))],[preIcClustIdx;interIcClustIdx]);

fprintf('For %s, there are\n %d pre-ictal and\n %d interictal sequences\n\n\n',...
    pt(whichPt).name,length(preIcClustIdx), length(interIcClustIdx));
fprintf(['For %s, regarding whether the pre-ictal period\n has a different cluster'...
    ' distribution from the interictal period,\n the p-value is %1.1e\n\n'],pt(whichPt).name,p_2);

%% Validate the chi2 with bootstrap
if 1 == 0
truePreIcIdx = preIcIdx; % the real pre ictal indices
trueInterIcIdx = interIcIdx; % the real interictal indices
truePreIcClustIdx = preIcClustIdx;
trueInterIcClustIdx = interIcClustIdx;

nPreIc = length(preIcIdx);
nInterIc = length(interIcIdx);
allClustIdx = [preIcClustIdx;interIcClustIdx];

% My goal is going to be to randomly swap out pre-ictal and interictal
% indices
nboot = 1e4;
chi2_boot = zeros(nboot,1);
for ib = 1:nboot
    p = randperm(length(allClustIdx),nPreIc);
    preIcClustIdx_temp = allClustIdx(p);
    interIcClustIdx_temp = allClustIdx(~ismember(1:end, p));
    [~,chi2_boot(ib)] = crosstab([ones(size(preIcClustIdx_temp));...
    2*ones(size(interIcClustIdx_temp))],...
    [preIcClustIdx_temp;interIcClustIdx_temp]);
end

[s_chi2_boot,I] = sort(chi2_boot);
diff_boot = abs(chi2_2-s_chi2_boot);
[~,I_min] = min(diff_boot);
boot_p = (length(I) - I_min)/length(I);

fprintf('By bootstrap, the p value is %1.1e\n',boot_p);
end

%% #2 Are 60 minute chunks containing seizures more likely to have certain cluster distributions
%[tbl_2,chi2_2,p_2,labels_2] = crosstab(sz_chunk,most_num);

%}

%% Save new pt struct
if allSpikes == 1
    save([saveFolder,'ptClustAll.mat'],'pt');
else
    save([saveFolder,'ptClust.mat'],'pt');
end

end




end
