function icAlone(pt,whichPts)

%{

I want to compare ranked data with binary data. I will have a
vector of all electrodes ranked by RL, where the earliest is #1 and those
that are not significantly involved are all tied for last. And then a
vector of all electrodes with 1 if in SOZ and 0 otherwise. This will check
to see if the earliest RL are more likely to be in SOZ than any other
electrodes, including both later RL as well as those electrodes that aren't
involved.

- One downside of this is that it may be getting a lot of its strength just
by looking at which electrodes spike the most. I would expect that the
electrodes that are spiking the most would be more likely to be within the
SOZ, if only because they are in the seizure...

Is there a way to combine the results of this for multiple
seizures/multiple patients?
- could I simply combine all the ranks for multiple seizures? So, for each
seizure, I would get the ranks of the SOZ electrodes and the ranks of the
non SOZ electrodes. I would lump all the SOZ ranks together across all
seizures and lump all the non-SOZ ranks together, and I would compare them.
There might be some independence issue but I can try it...

%}

%% Parameters
doClustOp = 0;
doLongPlots = 0;
doPlots = 1;

%% Ideal cluster numbers
ideal_k(1) = 3; %HUP064
ideal_k(2) = 4; %HUP065
ideal_k(3) = 3; %HUP068
ideal_k(6) = 3; %HUP074
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

%% Remove sequences with too many ties
if 1 == 0
keep = ones(size(seq_matrix,2),1);
old_seq_matrix = seq_matrix;
for s = 1:size(seq_matrix,2)
   curr_seq = seq_matrix(:,s); 
   curr_seq = curr_seq(~isnan(curr_seq));
   norepeats = unique(curr_seq);
   if length(norepeats) < 0.5*length(curr_seq)
       keep(s) = 0;
   end
end

seq_matrix(:,keep==0) = [];
fprintf('Deleted %d (%1.1f of all sequences) for containing >50 percent ties\n',...
    sum(keep==0),sum(keep==0)/length(keep));

end

%% Get lead channels and locations
[all_times,lead] = min(seq_matrix,[],1);
locs = pt(whichPt).electrodeData.locs(:,2:4);
lead_locs = locs(lead,:);

main_lead = mode(lead);

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

%% Break it down by seizure
inSOZ = [];
outSOZ = [];
for j =1:length(pt(whichPt).sz)
    
    seqs = pt(whichPt).sz(j).seq_matrix;
    [times,l] = min(seqs,[],1);
    sz_lead = mode(l);
    soz = pt(whichPt).sz(j).chs;
    if ismember(0,soz) == 1
        fprintf('warning, could not find all soz electrodes\n');
        soz(soz==0) = [];
    end
    
    % Discretize into 10 chunks to see how seizure changes
    nchunks = 1;
    [Y,~] = discretize(times,nchunks);
    
    figure
    % darker means more
    gs = (gray(50));
    for t = 1:nchunks
        
        % Get sequences in that chunk
        tseqs = seqs(:,find(Y==t));
        
        % Get average recruitment latency of each channel for that chunk
        tlat = nanmean(tseqs - min(tseqs,[],1),2);
        
        % get number of sequences in that chunk for each channel
        nseq = sum(~isnan(tseqs),2);
        
        
        % Find channels that are involved in sequences enough
        lambda = sum(sum(~isnan(tseqs)))/size(tseqs,1); 
        %all spikes divided by channels gives the expected number of spikes
        %per channel
        
        alpha1 = 95;
        X = poissinv(alpha1/100,lambda);
        sig_ch = find(nseq>X);
        
        % recruitment latency of significant channels
        tlat_sig = tlat(sig_ch);
        [Ylat,~] = discretize(tlat_sig,size(gs,1));
        [~,earliest] = min(tlat_sig);
        early_sig = sig_ch(earliest);
        
        % Not used
        if 1 == 0
        [~,most_common] = max(nseq); 
        [~,cl] = min(tseqs,[],1);
        chunk_lead = mode(cl);
        [Ycol,~] = discretize(nseq,size(gs,1));
        end
        
        scatter3(locs(:,1),locs(:,2),locs(:,3),60,'k')
        hold on
        scatter3(locs(sig_ch,1),locs(sig_ch,2),locs(sig_ch,3),...
            60,gs(Ylat,:),'filled');
        %scatter3(locs(Ylat(~isnan(Ylat)),1),locs(Ylat(~isnan(Ylat)),2),...
        %    locs(Ylat(~isnan(Ylat)),3),60,gs(Ylat(~isnan(Ylat)),:),'filled');
        %scatter3(locs(:,1),locs(:,2),locs(:,3),60,gs(Ycol,:),'filled');
        %scatter3(locs(nseq>X,1),locs(nseq>X,2),locs(nseq>X,3),60,'r','filled');
        %scatter3(locs(most_common,1),locs(most_common,2),locs(most_common,3),...
        %    60,'b','filled');
        
        scatter3(locs(soz,1),locs(soz,2),locs(soz,3),200,'r+');
        scatter3(locs(early_sig,1),locs(early_sig,2),locs(early_sig,3),...
            100,'g+');
        %{
        if isnan(chunk_lead) == 0
            scatter3(locs(chunk_lead,1),locs(chunk_lead,2),locs(chunk_lead,3),...
        100,'g+')
        end
        %}
        
        %% Check alignment of recruitment latency with SOZ
        tlat(nseq<X) = Inf;
        [~,I] = sort(tlat);
        r = 1:length(tlat);
        r(I) = r;
        r = r';
        
        % This is the rank of all electrodes by recruitment latency
        r(r>sum(nseq>=X)) = sum(nseq>=X) + 1;
        
        issoz = zeros(size(locs,1),1);
        issoz(soz) = 1;
        
        [p,h,stats] = ranksum(r(issoz==1),r(issoz==0));
        fprintf('P value is %1.1e for RL\n',p);
        
        %{
        inSOZ = [inSOZ;r(issoz==1)];
        outSOZ = [outSOZ;r(issoz==0)];
        %}
        
        
        fake_r = zeros(size(locs,1),1);
        fake_r(nseq>=X) = 1;
        %[p2,h2,stats2] = ranksum(fake_r(issoz==1),fake_r(issoz==0));
        %fprintf('P value is %1.1e for frequent spikes\n',p2);
        
        
        inSOZ = [inSOZ;fake_r(issoz==1)];
        outSOZ = [outSOZ;fake_r(issoz==0)];
        
        % Just look at the frequent spikers and re-run RL analysis
        if 1 == 0
        r_freq = r;
        issoz_freq = issoz;
        r_freq(nseq<X) = [];
        issoz_freq(nseq<X) = [];
        [p3,h3,stats3] = ranksum(r_freq(issoz_freq==1),r_freq(issoz_freq==0));
        fprintf('P value is %1.1e for frequent spikes\n',p3);
        end
        
        %% Do nearest SOZ analysis
        % Get lowest rlat of significant channels
        %minSOZ
        
        pause
        close(gcf)
       
    end
    
   
end
 
%%  combining all seizures
[p4,h4,stats4] = ranksum(inSOZ,outSOZ);
fprintf(['When I lump together all seizures, the p value for whether\n'...
    'SOZ electrodes have a different RL rank from non-SOZ electrodes is\n'...
    'p = %1.1e.\n'],p4);

%% Investigate spike timing

% Cluster spikes together
t_clust = all_times(1);
for i = 2:length(all_times)
    if all_times(i) - t_clust(end) < 0.3
        continue;
    else
        t_clust(end+1) = all_times(i);
    end
    
end
time_diff = diff(t_clust);
time_diff = [0 time_diff];

%plot(time_diff);
if 1 == 0
figure
scatter(t_clust,time_diff)
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
    if 1 == 1
    for k = 1:size(C,1)
    scatter3(C(k,1),C(k,2),C(k,3),60,colors(k,:),'filled');
    end
    end
    if 1 == 0
    scatter3(locs(main_lead,1),locs(main_lead,2),locs(main_lead,3),...
        60,'r','filled')
    end
    scatter3(locs(pt(whichPt).newSOZChs,1),locs(pt(whichPt).newSOZChs,2),...
        locs(pt(whichPt).newSOZChs,3),20,'k','filled')
     pause
     close(gcf)
end
    



%save([saveFolder,'ptClust.mat'],'pt');
    
end




end