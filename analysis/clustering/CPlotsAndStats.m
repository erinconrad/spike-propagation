
function CPlotsAndStats(pt,whichPts)

%% Parameters
window = 3600;
nboot = 1e4;

% Save file location
[~,~,~,resultsFolder,~] = fileLocations;

for whichPt = whichPts

fprintf('Doing %s\n',pt(whichPt).name);
if isempty(pt(whichPt).electrodeData) == 1
    continue
end

xyChan = pt(whichPt).electrodeData.locs;
saveFolder = [resultsFolder,'plots/',pt(whichPt).name,'/','clusters/'];
mkdir(saveFolder)


    
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
all_times(bad_idx) = [];
all_seq_cat(:,bad_idx) = [];
cluster_vec(bad_idx,:) = [];
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

%% Plot the clustering algorithm results 
figure
scatter3(cluster_vec(:,1),cluster_vec(:,2),cluster_vec(:,3),60,c_idx)
hold on
for i = 1:size(C,1)
   if i == bad_cluster, continue; end
   scatter3(C(i,1),C(i,2),C(i,3),100,colors(i,:),'filled');
     
end

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
figure
set(gcf,'Position',[50 100 1200 700])

% Sequence frequency
subplot(3,1,1)
[t_return,counts] = binCounts(all_times,window);
plot(t_return/3600,counts,'k','LineWidth',3);
hold on
for j = 1:length(pt(whichPt).sz)
   yl = ylim; 
   szOnset = pt(whichPt).sz(j).onset;
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
scatter(all_times/3600,cluster_vec(:,i)+repmat(toAdd,size(cluster_vec,1),1),20,c_idx)
hold on
ytick_locations(i) = toAdd+median(cluster_vec(:,i));
%text(all_times(1)/3600-0.3,toAdd+median(cluster_vec(:,i)),sprintf('%s',ttext{i}),'FontSize',30);
if i ~=3
    toAdd = toAdd + 20+(max(cluster_vec(:,i)) - min(cluster_vec(:,i+1)));%quantile(firstChs(:,i),0.95) - quantile(firstChs(:,i+1),0.05);
end
end
yl = ylim;
ylim([min(cluster_vec(:,1)), yl(2)])
for j = 1:length(pt(whichPt).sz)
   yl = ylim; 
   szOnset = pt(whichPt).sz(j).onset;
   sz = plot([szOnset szOnset]/3600,yl,'k','LineWidth',3);
end
yticks(ytick_locations)
yticklabels({'X','Y','Z'})
xlim([all_times(1)/3600-1 all_times(end)/3600+1])
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
    scatter(all_times/3600,cluster_vec(:,i)+repmat(toAdd,size(cluster_vec,1),1),20,c_idx)
    hold on
    ytick_locations(i-3) = toAdd;
    %text(all_times(1)/3600-0.3,toAdd,sprintf('%s',ttext{i-3}),'FontSize',30);
    if i ~=6
        toAdd = toAdd + 3;%quantile(final_vecs(:,i),0.9999) - quantile(final_vecs(:,i+1),0.0001); 
    end
end

yticks(ytick_locations)
yticklabels({'X','Y','Z'})
for j = 1:length(pt(whichPt).sz)
   yl = ylim; 
   szOnset = pt(whichPt).sz(j).onset;
   sz = plot([szOnset szOnset]/3600,yl,'k--','LineWidth',3);
end
xlim([all_times(1)/3600-1 all_times(end)/3600+1])
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
clust{i} = all_times(idx == i);
end

[sum_c,sum_times] = movingSumCounts(clust,all_times,window);

totalSum = zeros(1,size(sum_times,2));
for i = clusters
    totalSum = totalSum + sum_c(i,:);
end

subplot(3,1,3)

pl = zeros(k,1);
for i = clusters
   pl(i)= plot(sum_times/3600,sum_c(i,:)./totalSum,'color',colors(i,:),'LineWidth',3);
%plot(chunk_times,n_clusters_chunk(:,i),'color',colors(i,:),'LineWidth',2);
hold on
end

for j = 1:length(pt(whichPt).sz)
   yl = ylim; 
   szOnset = pt(whichPt).sz(j).onset;
   sz = plot([szOnset szOnset]/3600,yl,'k','LineWidth',3);
end
xlim([all_times(1)/3600-1 all_times(end)/3600+1])

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
destFolder = [resultsFolder,'pretty_plots/Fig2/'];
print(gcf,[destFolder,'clustTime'],'-depsc');
eps2pdf([destFolder,'clustTime.eps'])


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
num_cluster = zeros(n_chunks,k);
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


%% Validate the above method with bootstrap
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

%% Does sequence frequency predict sz?
time1 = sum(interIcTimes(:,2) - interIcTimes(:,1));
time2 = sum(preIcTime(:,2) - preIcTime(:,1));
count1 = length(interIcClustIdx);
count2 = length(preIcClustIdx);
[p_freq,z_freq] = poisson_mean_diff(count1,count2,time1,time2);

fprintf(['For %s, regarding whether the pre-ictal period\n has a different sequence'...
    ' frequency from the interictal period,\n the p-value is %1.1e\n\n'],...
    pt(whichPt).name,p_freq);


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

