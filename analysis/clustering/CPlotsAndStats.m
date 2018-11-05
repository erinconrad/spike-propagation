
function CPlotsAndStats(pt,whichPts)

%% Parameters
window = 3600;

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


%% Assign each sequence a color based on what cluster index it it
colors = [0 0 1;1 0 0;0 1 0; 0.5 0.5 1; 1 0.5 0.5; 0.5 1 0.5; 0.4 0.7 0.4];
c_idx = zeros(size(idx,1),3);
for i = 1:length(idx)
   c_idx(i,:) = colors(idx(i),:); 
end

if 1 == 1
%% Do plot
figure
set(gcf,'Position',[50 100 1200 1200])

% Plot of x, y, z coordinates of starting position over time
subplot(4,1,1)
toAdd = 0;
%marker = {'x','o','>'};
ttext = {'x','y','z'};
for i = 1:3
scatter(all_times/3600,cluster_vec(:,i)+repmat(toAdd,size(cluster_vec,1),1),20,c_idx)
hold on
text(all_times(1)/3600-0.3,toAdd+median(cluster_vec(:,i)),sprintf('%s',ttext{i}),'FontSize',30);
if i ~=3
    toAdd = toAdd + 10+(max(cluster_vec(:,i)) - min(cluster_vec(:,i+1)));%quantile(firstChs(:,i),0.95) - quantile(firstChs(:,i+1),0.05);
end
end
for j = 1:length(pt(whichPt).sz)
   yl = ylim; 
   szOnset = pt(whichPt).sz(j).onset;
   sz = plot([szOnset szOnset]/3600,yl,'k','LineWidth',2);
end
set(gca,'ytick',[]);
xlim([all_times(1)/3600-1 all_times(end)/3600+1])
title(sprintf('X, y, z coordinates of spike leader for %s',pt(whichPt).name));
set(gca,'FontSize',15);

% Plot of x, y, z coordinates of unit vector over time
subplot(4,1,2)
toAdd = 0;
%marker = {'x','o','>'};
for i = 4:6
    scatter(all_times/3600,cluster_vec(:,i)+repmat(toAdd,size(cluster_vec,1),1),20,c_idx)
    hold on
    
    text(all_times(1)/3600-0.3,toAdd,sprintf('%s',ttext{i-3}),'FontSize',30);
    if i ~=6
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

subplot(4,1,3)

pl = zeros(k,1);
for i = clusters
   pl(i)= plot(sum_times/3600,sum_c(i,:)./totalSum,'color',colors(i,:),'LineWidth',2);
%plot(chunk_times,n_clusters_chunk(:,i),'color',colors(i,:),'LineWidth',2);
hold on
end

for j = 1:length(pt(whichPt).sz)
   yl = ylim; 
   szOnset = pt(whichPt).sz(j).onset;
   sz = plot([szOnset szOnset]/3600,yl,'k','LineWidth',2);
end
xlim([all_times(1)/3600-1 all_times(end)/3600+1])

leg_text = {};
for i = clusters
   leg_text = [leg_text,sprintf('Cluster %d',i)];
end

legend([pl(pl~=0)],...
    leg_text,'Position',...
    [0.87 0.9 0.1 0.05]);
title(sprintf(['Proportion of sequences in given cluster, moving'...
    ' average %d s, %s'],...
    window,pt(whichPt).name));
set(gca,'FontSize',15);

%}


% Plot locations of centroids 


subplot(4,1,4)
scatter3(xyChan(:,2),xyChan(:,3),xyChan(:,4),60,'k');
hold on

for k = 1:size(C,1)
    if ismember(k,bad_cluster), continue; end;
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
saveas(gcf,[saveFolder,pt(whichPt).name,'cluster.png']);
%close(gcf)

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

%% #2 Are 60 minute chunks containing seizures more likely to have certain cluster distributions
%[tbl_2,chi2_2,p_2,labels_2] = crosstab(sz_chunk,most_num);
end

end


end

