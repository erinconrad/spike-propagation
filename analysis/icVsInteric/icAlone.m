function icAlone(pt,whichPts)

%% Parameters
doClustOp = 0;
doLongPlots = 1;
doPlots = 1;

%% Ideal cluster numbers
ideal_k(3) = 3; %HUP068
ideal_k(8) = 4; %HUP078

% Save file location
[~,~,~,resultsFolder,~] = fileLocations;


for whichPt = whichPts
    
saveFolder = [resultsFolder,'ic/cluster_validation/',pt(whichPt).name,'/'];
mkdir(saveFolder)
    


%% Get sequences
seq_matrix = [];

for j = 1:length(pt(whichPt).sz)
    seq_matrix = [seq_matrix,pt(whichPt).sz(j).seq_matrix];    
end

%% Get lead channels and locations
[all_times,lead] = min(seq_matrix,[],1);
locs = pt(whichPt).electrodeData.locs(:,2:4);
lead_locs = locs(lead,:);



%% Elbow method
if doClustOp == 1
SSE = zeros(10,1);
% Test k 1-10
for k = 1:10
    SSE_temp = zeros(30,1);
    
    % Pick best of 30 runs
    for j = 1:30
        [idx_test,C_test] = kmeans(lead_locs,k);
        
        % Get SSE summing over all the clusters
        for i = 1:k
            SSE_temp(j) = SSE_temp(j) + sum(sum((lead_locs(idx_test == i,:) - ...
                repmat(C_test(i,:),size(lead_locs(idx_test == i,:),1),1)).^2));
            
        end
    end
    SSE(k) = min(SSE_temp);
end

figure
plot(1:10,SSE)
end

%% Do clustering algorithm
metric = zeros(30,1);
for i = 1:30
    [idx_all{i},C_all{i},sumd_all{i},D_all{i}] = ...
        kmeans(lead_locs,ideal_k(whichPt));
    metric(i) = sum(sumd_all{i});
end
% Take the results of the clustering algorithm that worked the best
[~,minidx] = min(metric);
idx = idx_all{minidx};
C = C_all{minidx};
D = D_all{minidx};

pt(whichPt).cluster.idx = idx;
pt(whichPt).cluster.C = C;
pt(whichPt).cluster.D = D;

%% Get and plot representative sequences

for i = 1:ideal_k(whichPt)
    [~,I] = sort(D(:,i));
    rep_seq{i} = seq_matrix(:,I(1:12));
end

pt(whichPt).cluster.rep_seq = rep_seq;
if doLongPlots == 1
    for i = 1:ideal_k(whichPt)
        outputFile = sprintf('seqs_cluster_%d',i);
        icShowSpecificSequences(pt,whichPt,rep_seq{i},1,outputFile)
    end
end

%% Do plots

% Get sz times
szTimes = pt(whichPt).newSzTimes;

% Prep all times for plotting
all_times_plot = all_times;
for i = 1:length(all_times)
   [~,closestSzIdx] = min(abs(all_times(i)-szTimes(:,1)));
   closestSzTime = szTimes(closestSzIdx,1);
   all_times_plot(i) = all_times(i) - closestSzTime+closestSzIdx*200;
    
end

% Assign each sequence a color based on what cluster index it it
colors = [0 0 1;1 0 0;0 1 0; 0.5 0.5 1; 1 0.5 0.5; 0.5 1 0.5; 0.4 0.7 0.4];
c_idx = zeros(size(idx,1),3);
for i = 1:length(idx)
   c_idx(i,:) = colors(idx(i),:); 
end

if doPlots == 1
   figure
   set(gcf,'Position',[50 100 1200 1200])
   subplot(2,1,1)
   % Plot of x, y, z coordinates of starting position over time
    
    toAdd = 0;
    ttext = {'x','y','z'};
    for i = 1:3
    scatter(all_times_plot,lead_locs(:,i)+repmat(toAdd,size(lead_locs,1),1),20,c_idx)
    hold on
    %text(all_times_plot(1)-0.3,toAdd+median(lead_locs(:,i)),sprintf('%s',ttext{i}),'FontSize',30);
    if i ~=3
        toAdd = toAdd + 50+(max(lead_locs(:,i)) - min(lead_locs(:,i+1)));%quantile(firstChs(:,i),0.95) - quantile(firstChs(:,i+1),0.05);
    end
    end
    
    subplot(2,1,2)
    scatter3(locs(:,1),locs(:,2),locs(:,3),60,'k');
    hold on
    for k = 1:size(C,1)
    scatter3(C(k,1),C(k,2),C(k,3),60,colors(k,:),'filled');
    end
end
    



save([saveFolder,'ptClust.mat'],'pt');
    
end




end