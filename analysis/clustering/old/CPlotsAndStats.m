
function pt = CPlotsAndStats(pt,whichPts)

%% Parameters
allSpikes =0;
window = 3600;
nboot = 1e4;

% Save file location
[~,~,~,resultsFolder,~] = fileLocations;
destFolder = [resultsFolder,'clustering/'];
mkdir(destFolder);

allCounts = [];
allPat = [];
allChunk = [];
chi_tables = cell(max(whichPts),1);
chi_tables2 = cell(max(whichPts),1);
chi_tables_plot = cell(max(whichPts),1);

for whichPt = whichPts
    

fprintf('Doing %s\n',pt(whichPt).name);
if isempty(pt(whichPt).electrodeData) == 1
    continue
end

xyChan = pt(whichPt).electrodeData.locs;
saveFolder = [resultsFolder,'plots/',pt(whichPt).name,'/','clusters/'];
mkdir(saveFolder)


%% Get sz times and SOZ
locs = xyChan(:,2:4);
szTimes = pt(whichPt).newSzTimes;
soz = pt(whichPt).newSOZChs;
    
%% Pull cluster info
all_times = pt(whichPt).cluster.all_times;
all_seq_cat = pt(whichPt).cluster.all_seq_cat;
cluster_vec = pt(whichPt).cluster.cluster_vec;
k = pt(whichPt).cluster.k;
idx = pt(whichPt).cluster.idx;
C = pt(whichPt).cluster.C;
bad_cluster = pt(whichPt).cluster.bad_cluster;


%% Remove bad clusters
bad_idx = find(ismember(idx,bad_cluster));
if allSpikes == 1
    all_times_all = pt(whichPt).cluster.all_times_all;
    all_spikes = pt(whichPt).cluster.all_spikes;
    all_locs = pt(whichPt).cluster.all_locs;
    all_times_all(bad_idx) = [];
    all_locs(bad_idx,:) = [];
    all_spikes(bad_idx) = [];
    
else
    all_times(bad_idx) = [];
    all_seq_cat(:,bad_idx) = [];
    cluster_vec(bad_idx,:) = [];
end
idx(bad_idx) = [];
clusters = 1:k; clusters(bad_cluster) = [];

%% Remove extra times
if size(all_times,2) > size(cluster_vec,1)
    fprintf('Warning, there are %d more times than vectors\n',size(all_times,2)-size(cluster_vec,1));
    all_times = all_times(1:end-(size(all_times,2)-size(cluster_vec,1)));
end

%% Assign each sequence a color based on what cluster index it it
colors = [0 0 1;1 0 0;0 1 0; 0.5 0.5 1; 1 0.5 0.5; 0.5 1 0.5; 0.4 0.7 0.4];
c_idx = zeros(size(idx,1),3);
for i = 1:length(idx)
   c_idx(i,:) = colors(idx(i),:); 
end

%% Get sz onset chs
%{
szChs = [];
szIndex = [];
for j = 1:length(pt(whichPt).sz)
    szChs = [szChs;pt(whichPt).sz(j).chs]; 
    szIndex = [szIndex;j*ones(length(pt(whichPt).sz(j).chs),1)];
end
szByCh = cell(max(szChs),1);
for i = 1:length(szChs)
   szByCh{szChs(i)} = [szByCh{szChs(i)} szIndex(i)];
end
%}

if 1 == 1
%% Do plots

%% Centroid locations compared to SOZ

figure
scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
hold on

for i = 1:size(C,1)
   if i == bad_cluster, continue; end
   scatter3(C(i,1),C(i,2),C(i,3),100,colors(i,:),'filled');
     
end
scatter3(locs(soz,1),locs(soz,2),locs(soz,3),40,'k','filled')

%{
figure
scatter3(cluster_vec(:,4),cluster_vec(:,5),cluster_vec(:,6),60,c_idx)
hold on
for i = 1:size(C,1)
   if i == bad_cluster, continue; end
   scatter3(C(i,4),C(i,5),C(i,6),100,colors(i,:),'filled');
     
end
%}

%% Main cluster plot
if allSpikes == 1
   plot_thing = all_locs;
   plot_times = all_times_all;
else
   plot_thing = cluster_vec;
   plot_times = all_times;
end

figure
set(gcf,'Position',[50 100 1200 700])

% Sequence frequency
subplot(3,1,1)
[t_return,counts] = binCounts(plot_times,window);
plot(t_return/3600,counts,'k','LineWidth',3);
hold on
for j = 1:size(szTimes,1)
   yl = ylim; 
   szOnset = szTimes(j,1);
   sz = plot([szOnset szOnset]/3600,yl,'k','LineWidth',3);
end
xlim([all_times(1)/3600-1 all_times(end)/3600+1])
ylabel('Sequences per hour');
title(sprintf('Sequence frequency'));
set(gca,'FontSize',20);

% Plot of x, y, z coordinates of starting position over time
subplot(3,1,2)
toAdd = 0;
%marker = {'x','o','>'};
ttext = {'x','y','z'};
ytick_locations = zeros(3,1);
for i = 1:3
scatter(plot_times/3600,plot_thing(:,i)+repmat(toAdd,size(plot_thing,1),1),20,c_idx)
hold on
ytick_locations(i) = toAdd+median(plot_thing(:,i));
%text(all_times(1)/3600-0.3,toAdd+median(cluster_vec(:,i)),sprintf('%s',ttext{i}),'FontSize',30);
if i ~=3
    toAdd = toAdd + 20+(max(plot_thing(:,i)) - min(plot_thing(:,i+1)));%quantile(firstChs(:,i),0.95) - quantile(firstChs(:,i+1),0.05);
end
end
yl = ylim;
ylim([min(plot_thing(:,1)), yl(2)])
for j = 1:size(szTimes,1)
   yl = ylim; 
   szOnset = szTimes(j,1);
   sz = plot([szOnset szOnset]/3600,yl,'k','LineWidth',3);
end
yticks(ytick_locations)
yticklabels({'X','Y','Z'})
xlim([plot_times(1)/3600-1 plot_times(end)/3600+1])
ylabel('Coordinate');
title(sprintf('X, Y, Z coordinates of spike leader'));
set(gca,'FontSize',20);

if 1 == 0
% Plot of x, y, z coordinates of unit vector over time
subplot(4,1,3)
toAdd = 0;
%marker = {'x','o','>'};
ytick_locations = zeros(3,1);
for i = 4:6
    scatter(plot_times/3600,cluster_vec(:,i)+repmat(toAdd,size(cluster_vec,1),1),20,c_idx)
    hold on
    ytick_locations(i-3) = toAdd;
    %text(all_times(1)/3600-0.3,toAdd,sprintf('%s',ttext{i-3}),'FontSize',30);
    if i ~=6
        toAdd = toAdd + 3;%quantile(final_vecs(:,i),0.9999) - quantile(final_vecs(:,i+1),0.0001); 
    end
end

yticks(ytick_locations)
yticklabels({'X','Y','Z'})
for j = 1:size(szTimes,1)
   yl = ylim; 
   szOnset = szTimes(j,1);
   sz = plot([szOnset szOnset]/3600,yl,'k','LineWidth',3);
end
xlim([plot_times(1)/3600-1 plot_times(end)/3600+1])
ylabel('Coordinate');
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
for i = clusters
clust{i} = plot_times(idx == i);
end

[sum_c,sum_times] = movingSumCounts(clust,plot_times,window);

totalSum = zeros(1,size(sum_times,2));
for i = clusters
    totalSum = totalSum + sum_c(i,:);
end

subplot(3,1,3)

pl = zeros(k,1);
for i = clusters
   %pl(i)= plot(sum_times/3600,sum_c(i,:),'color',colors(i,:),'LineWidth',3);

   pl(i)= plot(sum_times/3600,sum_c(i,:)./totalSum,'color',colors(i,:),'LineWidth',3);
%plot(chunk_times,n_clusters_chunk(:,i),'color',colors(i,:),'LineWidth',2);
hold on
end

for j = 1:size(szTimes,1)
   yl = ylim; 
   szOnset = szTimes(j,1);
   sz = plot([szOnset szOnset]/3600,yl,'k','LineWidth',3);
end
xlim([plot_times(1)/3600-1 plot_times(end)/3600+1])

leg_text = {};
for i = clusters
   leg_text = [leg_text,sprintf('Cluster %d',i)];
end

legend([pl(pl~=0);sz],...
    [leg_text,'Seizures'],'Position',...
    [0.87 0.9 0.1 0.05],'FontSize',20);
xlabel('Time (hours)');
title(sprintf(['Proportion of sequences in given cluster']));
set(gca,'FontSize',20);

%saveas(gcf,[saveFolder,pt(whichPt).name,'cluster.png']);
%close(gcf)
print(gcf,[destFolder,'clustTime_',sprintf('%s',pt(whichPt).name)],'-depsc');
eps2pdf([destFolder,'clustTime_',sprintf('%s',pt(whichPt).name),'.eps'])


%% Plot locations of centroids 
if 1 == 0

figure
scatter3(xyChan(:,2),xyChan(:,3),xyChan(:,4),60,'k');
hold on



for i = 1:length(szByCh)
   if isempty(szByCh{i}) == 0
      szText = sprintf('\n');
      for j = 1:length(szByCh{i})
          szText = [szText,sprintf('%d, ',szByCh{i}(j))];
          
      end
      szText = szText(1:end-2);
      scatter3(xyChan(i,2),xyChan(i,3),xyChan(i,4),100,'b');
      text(xyChan(i,2),xyChan(i,3),xyChan(i,4),szText,'FontSize',10)
      
   end
    
end

%scatter3(xyChan(szChs,2),xyChan(szChs,3),xyChan(szChs,4),100,colors(szIndex,:));

for k = 1:size(C,1)
    if ismember(k,bad_cluster), continue; end
    scatter3(C(k,1),C(k,2),C(k,3),60,colors(k,:),'filled');
    plot3([C(k,1) C(k,1) + 10*C(k,4)],...
        [C(k,2) C(k,2) + 10*C(k,5)],...
        [C(k,3) C(k,3) + 10*C(k,6)],'k','LineWidth',2)
end

title(sprintf('Spike leader and propagation vectors for %s',pt(whichPt).name));
set(gca,'FontSize',15);
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'zticklabel',[])


end
%}





%% Statistical tests

% Get all seizure times
szAll = szTimes(:,1);


%% #1 does cluster distribution change over 60 minute chunks?

% 60 minute chunks
test_t = 3600;

% Determine most popular cluster
pop_c = mode(idx);

% Fix this to account for multiple time groups!!!!!!!!!
%{
totalTime = 0;
chunk_times = [];
for i = 1:size(pt(whichPt).allTimes,1)
    totalTime = totalTime + pt(whichPt).allTimes(2) - pt(whichPt).allTimes(1);  
end
n_chunks = ceil(totalTime/test_t);
%}


% Divide run into 60 minute chunks
n_chunks = ceil((plot_times(end) - plot_times(1))/test_t);

prop_pop = zeros(n_chunks,1);
times_pop = zeros(n_chunks,1);
which_chunk = zeros(length(all_times),1);
num_cluster = zeros(n_chunks,k);
sz_chunk = zeros(n_chunks,1);
most_num = zeros(n_chunks,1);

for i = 1:n_chunks
   curr_times = [plot_times(1)+(i-1)*test_t, min(plot_times(1) + i*test_t,plot_times(end))];
   
   % get the indices of the sequences in that time chunk
   curr_seqs = find(plot_times >= curr_times(1) & plot_times <= curr_times(2));
   prop_pop(i) = sum(idx(curr_seqs) == pop_c)/length(curr_seqs);
   
   % define the time chunk for those sequences
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

pt(whichPt).cluster.changetime.p = p_1;

%% Same thing, but just 2 most popular clusters and first 24 chunks
%{
most_common = mode(idx);
second_common = mode(idx(idx~=most_common));
new_idx = idx(idx == most_common | idx == second_common);
new_which_chunk = which_chunk(idx == most_common | idx == second_common);
unique_chunks = unique(new_which_chunk);
allowable_chunks = unique_chunks(1:24);
new_idx = new_idx(ismember(new_which_chunk,allowable_chunks));
new_which_chunk = new_which_chunk(ismember(new_which_chunk,allowable_chunks));

[tbl_5,chi2_5,p_5,labels_5] = crosstab(new_which_chunk,new_idx);
fprintf(['For %s, regarding whether 60 minute chunks\n have different cluster'...
    ' distributions,\n when just taking first 2 clusters and '...
    'first 24 hours,\nthe p-value is %1.1e\n\n\n'],pt(whichPt).name,p_1);

chi_tables2{whichPt} = tbl_5;
%}

%% Validate the above method with permutation
%{
% Randomly shuffle which chunk each sequence is in and recheck
chi2_1_test = zeros(nboot,1);
for ib = 1:nboot
    fake_chunks = zeros(length(idx),1);
    n_remain = length(idx);
    for i = 1:n_chunks
        y = randsample(n_remain,sum(which_chunk == i));
        fake_chunks(y) = i;
        n_remain = n_remain - length(y);
    end
    [~,chi2_1_test(ib)] = crosstab(fake_chunks,idx);  
end
[s_chi2_1_test,I] = sort(chi2_1_test);
diff_boot_1 = abs(chi2_1-s_chi2_1_test);
[~,I_min] = min(diff_boot_1);
boot_p_1 = (length(I) - I_min)/length(I);

fprintf('By bootstrap, the p value is %1.1e\n',boot_p_1);
%}

%% Do test to compare cluster distributions between pre-ictal and inter-ictal data

% Get all the pre-ictal sequences
preIcRange = [-60*60,-1*60];
preIcClustIdx = [];
preIcIdx = [];
preIcTimes = zeros(length(pt(whichPt).sz),2);
for j = 1:length(szAll)%1:length(pt(whichPt).sz)
    %szTime = pt(whichPt).sz(j).onset;
    
    % Get the current seizure time
    szTime = szAll(j);
    
    % Get the range of pre-ictal times (1 to 60 minutes before the seizure)
    preIcTime = szTime + preIcRange;
    preIcTimes(j,:) = preIcTime;
    
    % The cluster indices of the sequences falling within that range
    preIcClustIdx = [preIcClustIdx;idx(plot_times >= preIcTime(1) & ...
        plot_times <= preIcTime(2))];
    
    % the sequence numbers falling within that range
    if allSpikes == 1
        preIcIdx = [preIcIdx;find(plot_times >= preIcTime(1) & ...
        plot_times <= preIcTime(2))];
    else
        preIcIdx = [preIcIdx;find(plot_times >= preIcTime(1) & ...
            plot_times <= preIcTime(2))'];
    end
    
end

szTimes = szAll;



% Get inter-ictal sequences
n_chunks = length(szAll);
interIcTimes = [];
for i = 1:n_chunks
    
    % potential late_time is 3 hours before the seizure (2 hours before
    % pre-ictal period)
   late_time = preIcTimes(i,1) - 3600*2;
   
   % potential early time is either the start of the run or an hour after
   % the last seizure
   if i ==1
       early_time = plot_times(1);
   else
       early_time = szTimes(i-1) + 3600*1;
   end
   
   if early_time > late_time
       continue
   end
   
   interIcTimes = [interIcTimes;early_time late_time];

end
 % add time after last seizure (excluding last 2 hours in case there was a
 % seizure right after we stop recording)
 late_time = plot_times(end)- 3600*2;
 early_time = szTimes(end) + 3600*1;
 if late_time>early_time
     interIcTimes = [interIcTimes;early_time late_time];
 end

interIcClustIdx = [];
interIcIdx = [];
for i = 1:size(interIcTimes,1)
    interIcClustIdx = [interIcClustIdx;idx(plot_times >= interIcTimes(i,1) & ...
        plot_times <= interIcTimes(i,2))];
    if allSpikes == 1
        interIcIdx = [interIcIdx;find(plot_times >= interIcTimes(i,1) & ...
        plot_times <= interIcTimes(i,2))];
    else
        
        interIcIdx = [interIcIdx;find(plot_times >= interIcTimes(i,1) & ...
            plot_times <= interIcTimes(i,2))'];
    end
    
end

[tbl_2,chi2_2,p_2,labels_2] = crosstab([ones(size(preIcClustIdx));...
    2*ones(size(interIcClustIdx))],[preIcClustIdx;interIcClustIdx]);

fprintf('For %s, there are\n %d pre-ictal and\n %d interictal sequences\n\n\n',...
    pt(whichPt).name,length(preIcClustIdx), length(interIcClustIdx));
fprintf(['For %s, regarding whether the pre-ictal period\n has a different cluster'...
    ' distribution from the interictal period,\n the p-value is %1.1e\n\n'],pt(whichPt).name,p_2);

chi_tables_plot{whichPt} = tbl_2;
pt(whichPt).cluster.pre_ic.p = p_2;

%% Get 2x2 table for just 2 most common clusters

%{
interIcClustIdx(interIcClustIdx~=most_common & interIcClustIdx~=second_common) = [];
preIcClustIdx(preIcClustIdx~=most_common & preIcClustIdx~=second_common) = [];

[tbl_3,chi2_3,p_3,labels_3] = crosstab([ones(size(preIcClustIdx));...
    2*ones(size(interIcClustIdx))],[preIcClustIdx;interIcClustIdx]);

fprintf('For %s, there are\n %d pre-ictal and\n %d interictal sequences\n\n\n',...
    pt(whichPt).name,length(preIcClustIdx), length(interIcClustIdx));
fprintf(['For %s, regarding whether the pre-ictal period\n has a different cluster'...
    ' distribution from the interictal period,\n'...
    'only taking into account 2 most common clusters,\nthe p-value is %1.1e\n\n'],pt(whichPt).name,p_3);

chi_tables{whichPt} = tbl_3;
%}

%% Does sequence frequency predict sz?
time1 = sum(interIcTimes(:,2) - interIcTimes(:,1));
time2 = sum(preIcTime(:,2) - preIcTime(:,1));
count1 = length(interIcClustIdx);
count2 = length(preIcClustIdx);
[p_freq,z_freq] = poisson_mean_diff(count1,count2,time1,time2);

fprintf(['For %s, regarding whether the pre-ictal period\n has a different sequence'...
    ' frequency from the interictal period,\n the p-value is %1.1e\n\n'],...
    pt(whichPt).name,p_freq);

allCounts = [allCounts;count1/time1;count2/time2];
allChunk = [allChunk;0;1];
allPat = [allPat;whichPt; whichPt];

%% Validate the chi2 with bootstrap
if 1==0
truePreIcIdx = preIcIdx; % the real pre ictal indices
trueInterIcIdx = interIcIdx; % the real interictal indices
truePreIcClustIdx = preIcClustIdx;
trueInterIcClustIdx = interIcClustIdx;

nPreIc = length(preIcIdx);
nInterIc = length(interIcIdx);
allClustIdx = [preIcClustIdx;interIcClustIdx];

% My goal is going to be to randomly swap out pre-ictal and interictal
% indices
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





%% Are the cluster leads close to the SOZ electrodes?
if 1 == 0

%% NEED TO DO SOMETHING TO DEAL WITH MULTIPLE MEASUREMENTS
% Like, if you get to pick the best of n randomly placed channels, how
% often will you get closer

szLoc = xyChan(unique(szChs),2:4);
szCentroid = mean(szLoc,1);

C_mode = C(pop_c,:);
D_mode = sqrt(sum((C_mode(1:3)-szCentroid).^2));

% Decide how to remove bad clusters!

C(bad_cluster,:) = [];

cToSz = zeros(size(C,1),size(szLoc,1));
for i = 1:size(szLoc,1)
   for j = 1:size(C,1) 
       cToSz(j,i) = sqrt(sum((C(j,1:3)-szLoc(i,:)).^2));
    
   end
end

%{
figure
imagesc(cToSz)
xlabel('Seizure onset zone electrode')
ylabel('Cluster centroid')
colorbar
%}


avgD = zeros(size(C,1),1);
for i = 1:size(C,1)
   %avgD(i) = mean(cToSz(i,:)); 
   avgD(i) = sqrt(sum((C(i,1:3)-szCentroid).^2));
end

allChsAvgD = zeros(size(xyChan,1),1);
allChToSz = zeros(size(xyChan,1),size(szLoc,1));
for i = 1:size(szLoc,1)
   for j = 1:size(xyChan,1) 
       allChToSz(j,i) = sqrt(sum((xyChan(j,2:4)-szLoc(i,:)).^2));
    
   end
end
%{
figure
imagesc(allChToSz)
xlabel('Seizure onset zone electrode')
ylabel('Channel locations')
colorbar
%}

for i = 1:size(xyChan,1)
   %allChsAvgD(i) = mean(allChToSz(i,:)); 
   allChsAvgD(i) = sqrt(sum((xyChan(i,2:4)-szCentroid).^2));
end

[bestD,bestC] = min(avgD);
[allChsS,I] = sort(allChsAvgD);
diff_boot = abs(bestD-allChsS);
[~,I_min] = min(diff_boot);
perc = (length(I) - I_min)/length(I);
fprintf(['The best cluster of %d is\n %1.1f mm from the average SOZ electrode location,\n'...
    'closer to the SOZ electrodes than \n'...
    '%1.2f (%d of %d) of all electrodes\n\n\n'],...
    size(C,1),bestD,perc,length(I) - I_min,length(I));


diff_boot = abs(D_mode-allChsS);
[~,I_min] = min(diff_boot);
perc = (length(I) - I_min)/length(I);
fprintf(['The mode cluster is\n %1.1f mm from the average SOZ electrode location,\n'...
    'closer to the SOZ electrodes than \n'...
    '%1.2f (%d of %d) of all electrodes\n\n\n'],D_mode,perc,length(I) - I_min,length(I));


end
end


end

%% Bar graphs
figure
[ha, pos] = tight_subplot(1,length(whichPts),[.01 .01],[.1 .08],[.02 .01]); 
for j = 1:length(whichPts)
    axes(ha(j));
    tbl = chi_tables_plot{whichPts(j)};
    new_tbl = tbl;
    for k = 1:size(new_tbl,1)
       new_tbl(k,:) = new_tbl(k,:)/sum(new_tbl(k,:));
        
    end
    bar(new_tbl)
    xticklabels({'Pre-ic','Inter-ic'});
    
    legend_names = cell(size(tbl,2),1);
    for k = 1:length(legend_names)
       legend_names{k} = sprintf('Cluster %d',k);  
    end
    yticklabels([])
    if j == 4
        lgnd=legend(legend_names,'location','northwest');
        set(lgnd,'color','none');
    end
    title(sprintf('%s',pt(whichPts(j)).name));
    if j == 1
        ylabel('Proportion of sequences');
    end
    set(gca,'FontSize',23)
end
%pause
fig = gcf;
fig.PaperUnits = 'inches';
posnow = get(fig,'Position');
print(gcf,[destFolder,'chi2_allspikes'],'-dpng');

%{
%% Get full table for chi_2 to put into R
names = cell(length(whichPts)*4,1);
cluster = cell(length(whichPts)*4,1);
state = cell(length(whichPts)*4,1);
counts = zeros(length(whichPts)*4,1);

count = 1;
for whichPt = whichPts
    names{count} = pt(whichPt).name;
    names{count+1} = pt(whichPt).name;
    names{count+2} = pt(whichPt).name;
    names{count+3} = pt(whichPt).name;
    
    cluster{count} = 'Cluster1';
    cluster{count+1} = 'Cluster2';
    cluster{count+2} = 'Cluster1';
    cluster{count+3} = 'Cluster2';
    
    state{count} = 'Preic';
    state{count+1} = 'Preic';
    state{count+2} = 'Interic';
    state{count+3} = 'Interic';
    
    
    % Get appropriate table
    tbl = chi_tables{whichPt};
    clust1_preic = tbl(1,1);
    clust2_preic = tbl(1,2);
    clust1_interic = tbl(2,1);
    clust2_interic = tbl(2,2);
    
    counts(count) = clust1_preic;
    counts(count+1) = clust2_preic;
    counts(count+2) = clust1_interic;
    counts(count+3) = clust2_interic;
    
    count = count + 4;
end

T = table(names,state,cluster,counts);



%% Get full table for chi_2 to put into R for cluster dist over time
names = cell(length(whichPts)*48,1);
cluster = cell(length(whichPts)*48,1);
state = zeros(length(whichPts)*48,1);
counts = zeros(length(whichPts)*48,1);


count = 1;
for whichPt = whichPts
    for j = 0:47
       names{count+j} = pt(whichPt).name; 
    end
    
    for j = 0:2:46
       cluster{count+j} = 'Cluster1'; 
    end
    for j = 1:2:47
       cluster{count+j} = 'Cluster2'; 
    end
    
    for j = 0:47
       state(count+j) = floor((1+j+1)/2);
    end
    
    tbl = chi_tables2{whichPt};
    tbl = tbl';
    tbl_single = tbl(:);
    
    for j = 0:47
       counts(count+j) = tbl_single(j+1);
    end
    
    count = count+48;
    
end

T = table(names,state,cluster,counts);
%}


end

