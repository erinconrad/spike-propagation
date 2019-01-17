function CInfluence(pt,cluster,whichPts)

%{

SUPER important to make sure we align the spikes here!!!!!!
SHould double check

%}

% Parameters
doPlots = 0; %0 = no, 1=normal, 2=pretty
plotConn = 0;
removeTies = 1;
doBootstrap = 0;
alpha1 = 95;
map_text = 'jet';
fh_map = str2func(map_text);

[~,~,scriptFolder,resultsFolder,~] = fileLocations;
p1 = genpath(scriptFolder);
addpath(p1);

destFolder = [resultsFolder,'influence/'];
mkdir(destFolder)

% Initialize matrices for all patients
allFreqDist = [];
allAllDist =[];
allSADist = [];
allSpikeDist = [];
allSAToFreqDist = [];
allSAs = [];
allSpikers = [];
allSF = [];
allDegPref = [];
allMaxSAs = [];

if isempty(whichPts) == 1
    for i = 1:length(pt)
        if isempty(pt(i).seq_matrix) == 0
            if size(cluster(i).bad_cluster) < cluster(i).k
                whichPts = [whichPts,i];
            end
        end
    end
elseif whichPts == 100
    whichPts = [1,4,6,7,8,9,12,14,15,16,17,18,19,20,22,24,25,27,30,31];
elseif whichPts == 300
    whichPts = [1 4 6 8 9 12 17 18 19 20 22 24 25 27 30 31];
end

if isequal(whichPts,[1,4,6,7,8,9,12,14,15,16,17,18,19,20,22,24,25,27,30,31]) == 0
    error('Warning, not doing correct patients!\n');
end

for whichPt = whichPts
    
    %% Patient parameters
    fprintf('Doing %s\n',pt(whichPt).name);
    locs = pt(whichPt).electrodeData.locs(:,2:4);
    nchs = size(locs,1);
    szTimes = pt(whichPt).newSzTimes;
    soz = pt(whichPt).newSOZChs; 
    saveFolder = [destFolder,pt(whichPt).name,'/'];
    mkdir(saveFolder);
    
    if isempty(soz) == 1
        fprintf('WARNING, soz empty for %s, skipping\n',pt(whichPt).name);
        continue
    end
    
    seq_matrix = pt(whichPt).seq_matrix;
    
    %% Remove ties
    if removeTies == 1
        keep = ones(size(seq_matrix,2),1);
        for s = 1:size(seq_matrix,2)
           curr_seq = seq_matrix(:,s);
           nonans = curr_seq(~isnan(curr_seq));
           norepeats = unique(nonans);
           if length(norepeats) < 0.5*length(nonans)
               keep(s) = 0;
           end
        end
        seq_matrix(:,keep==0) = [];
        fprintf(['%s had %d sequences (%1.2f of all sequences) deleted'...
        'for having >50 percent ties\n%d sequences remain\n'],...
        pt(whichPt).name,sum(keep == 0),sum(keep == 0)/length(keep),sum(keep==1));

    end
    
    %% Remove ictal sequences
    all_times = seq_matrix(:);
    icTimes = find(any(all_times >= (szTimes(:,1)-repmat(60,size(szTimes,1),1))' ...
        & all_times <= szTimes(:,2)',2));
    seq_matrix(icTimes) = nan;
    fprintf('Removed %d ictal spikes\n',length(icTimes));
    %{
    first_time = min(seq_matrix,[],1);
    t = (any(first_time >= (szTimes(:,1)-repmat(60,size(szTimes,1),1)) ...
        & first_time <= szTimes(:,2),2));
    seq_matrix(:,t) = [];
    fprintf('Removed %d ictal spikes \n',sum(t));
    %}
    
    %% Get cluster info
    all_times_all = cluster(whichPt).all_times_all; % all spike times
    all_spikes = cluster(whichPt).all_spikes; % all spike channels
    all_locs = cluster(whichPt).all_locs;
    k = cluster(whichPt).k; % the number of clusters
    idx = cluster(whichPt).idx; % the cluster index for every spike
    C = cluster(whichPt).C; % the centroids of the clusters
    bad_cluster = cluster(whichPt).bad_cluster; % which clusters are bad
    
    %% Compare number of spikes in cluster array and my data
    if sum(sum(~isnan(seq_matrix))) ~= length(all_times_all)
        error('Warning, number of spikes do not align\n');
    end
    
    
    
    %% Find bad spikes
    bad_idx = find(ismember(idx,bad_cluster));
    
    % Nx2 array of bad spikes, showing the channel and time
    bad_spikes = [all_spikes(bad_idx),all_times_all(bad_idx)];
    
    %% Get all sequences
    
    new_seq_matrix = seq_matrix;
    n_removed = 0;
    
    %% Go through sequence matrix and remove bad spikes
    for ich = 1:size(seq_matrix,1)
        % loop across electrodes
        
        % All spike times for this channel
        spikeTimesCh = seq_matrix(ich,:);
        
        % Get the bad spikes in that channel
        bad_times_for_ch = bad_spikes(bad_spikes(:,1) == ich,2);
        
        % Make sure I am finding all of them
        Lia = ismember(spikeTimesCh,bad_times_for_ch);
        if sum(Lia) ~= length(bad_times_for_ch)
            error(sprintf('Did not find all bad spikes for channel %d\n',ich));
        end
        
        %{
        if sum(Lia) > 0
            fprintf('Removed %d spikes for channel %d\n',sum(Lia),ich)
        end
        %}
        n_removed = n_removed + sum(Lia);
        
        % Make bad spikes nans
        spikeTimesCh(Lia==1) = nan;
        new_seq_matrix(ich,:) = spikeTimesCh;
        
        
    end
    
    if n_removed~=length(bad_idx)
        error('Incorrect number of bad spikes removed\n');
    end
    fprintf('Removed %d spikes for being in bad clusters\n',n_removed);
    
    
    %% Remove sequences that have fewer than 5 spikes
    removeSeq = zeros(size(new_seq_matrix,2),1);
    for s = 1:size(new_seq_matrix,2)
        currSeq = new_seq_matrix(:,s);
        currSeq(isnan(currSeq)) = [];
        if length(currSeq) < 5
            removeSeq(s) = 1;
        end
    end
    
    fprintf('Removed %d sequences for now being too short\n',sum(removeSeq));
    new_seq_matrix(:,removeSeq==1) = [];
    
    
    seq_matrix = new_seq_matrix;
    
    fprintf('%d sequences remain\n',size(seq_matrix,2));
    
    
    %% Construct a matrix of channel connections
    chCh = zeros(nchs,nchs);
    all_spike_times = [];
    upDown = zeros(nchs,2);
    for i = 1:size(seq_matrix,2)
        seq = seq_matrix(:,i);
        spike_chs = find(isnan(seq) == 0);
        spike_times = seq((isnan(seq) == 0));
        
        % sort channels by spike time
        [spike_times,I] = sort(spike_times);
        spike_chs = spike_chs(I);
        
        all_spike_times = [all_spike_times;spike_times];
        
        % Get all unidirectional pairwise connections that start with first
        % channel
        B = nchoosek(spike_chs,2);
        B = B(B(:,1) == spike_chs(1),:);
        
        for j = 1:size(B,1)
            chCh(B(j,1),B(j,2)) = chCh(B(j,1),B(j,2)) + 1;
        end
        
        % Add info for degree preference
        for j = 1:length(spike_chs)
            % Add the number of channels before it in the sequence
            upDown(spike_chs(j),1) = upDown(spike_chs(j),1) + j-1;
            
           upDown(spike_chs(j),2) = upDown(spike_chs(j),2) + length(spike_chs)-j;
        end
        
    end
    
    degPref = zeros(nchs,1);
    for i = 1:size(upDown,1)
        degPref(i) = 100*(upDown(i,2)-upDown(i,1))/(upDown(i,2)+upDown(i,1));
    end
    allDegPref = [allDegPref;degPref];
    
    if plotConn == 1
        figure
        imagesc(chCh)
        colorbar
    end
    
    if whichPt == 8, exampleChCh = chCh; end
    
    %% Do a test of symmetry
    % If spikes are randomly oriented, then chCh should be symmetric, which
    % is to say that for all i and j, chCh(i,j) ~= chCh(j,i)
    
    
    
    %% Get significant connections
    if doBootstrap == 1
        
        % Here, for each permutation, I am constructing a chCh matrix where
        % I am distributing the true total number of connections randomly
        % across all elements of the nch by nch matrix.
        ncons = sum(sum(chCh));
        nboot = 1e3;
        max_size = nchs*nchs;
        chCh_all = zeros(nboot,nchs,nchs);
        chCh_diff_all = zeros(nboot,nchs,nchs);
        for ib = 1:nboot
            if mod(ib,100) == 0
                fprintf('Doing %d of %d\n', ib,nboot);
            end
            chCh_f = zeros(nchs,nchs);
            for j = 1:ncons
               chCh_f(randi(max_size)) = chCh_f(randi(max_size)) + 1;
            end
            chCh_all(ib,:,:) = chCh_f;
            
            % Also calculate the difference between i,j and j,i
            for i = 1:size(chCh_f,1)
                for j = 1:size(chCh_f,2)
                    chCh_diff(i,j) = chCh_f(i,j) - chCh_f(j,i);
                end
            end
            chCh_diff_all(ib,:,:) = chCh_diff;
            
        end
        
        if 1 == 1
            figure
            imagesc(squeeze(mean(chCh_diff_all,1)))
            colorbar
        end
        
        % Get the 95% of number of connections
        s_con = sort(chCh_all(:));
        perc = prctile(s_con,alpha1);
       
           
        fprintf(['By permutation testing, the minimum number of counts for a\n'...
        'connection to be significant is\n'...
        '%d for an alpha of %1.1fth percentile\n\n'],perc,alpha1);
    
        % Get the 95% of differences 
        s_diff = sort(chCh_diff_all(:));
        perc_diff = prctile(s_diff,alpha1);
        
    end
    
    % Assume poisson distribution
    ncons = sum(sum(chCh));
    lambda = ncons/nchs^2;
    X = poissinv(alpha1/100,lambda);
    minCount = X;
    fprintf(['By poisson assumption, the number of counts is:\n'...
        '%d\n\n'],X);
    
    if whichPt == 8, exampleX = X; end
    
    %% Now find channels with more spikes than expected by chance.
    n_spikes = length(all_spike_times);
    lambda_spikes = n_spikes/nchs;
    X_spikes = poissinv(alpha1/100,lambda_spikes);
    minCountSpikes = X_spikes;
    
    n_spikes_ch = sum(~isnan(seq_matrix),2);
    ch_w_spikes = find(n_spikes_ch>minCountSpikes);
    spiker = n_spikes_ch>minCountSpikes;
    allSpikers = [allSpikers;spiker];
    
    %% Get the indices of the channels that are influenced by each other
    chInfluence = cell(nchs,1);
    for i = 1:length(chInfluence)
        for j = 1:size(chCh,2)
            if chCh(i,j) > minCount
                chInfluence{i} = [chInfluence{i},j];
            end
        end
    end
    
    %% Get the surface are of influence of each channel
    sa = zeros(nchs,1);
    for i = 1:nchs
        if isempty(chInfluence{i}) == 1, continue; end
        downstream = chInfluence{i};
        if length(downstream) <= 2, continue; end
        down_locs = locs(downstream,:);
        down_locs = [locs(i,:);down_locs];
        tri = delaunay(down_locs(:,1),down_locs(:,2));
        P = [down_locs(:,1),down_locs(:,2),down_locs(:,3)];
        v1 = P(tri(:,2),:)-P(tri(:,1),:);
        v2 = P(tri(:,3),:)-P(tri(:,2),:);
        cp = 0.5*cross(v1,v2);
        sa(i) = sum(sqrt(dot(cp,cp,2)));
    end
    
    
    
    
    
    %% How far is the channel with max area of influence from SOZ?
    allSF = [allSF;n_spikes_ch];
    allSAs = [allSAs;sa];
    allMaxSAs = [allMaxSAs;(max(sa))];
    
    % Distance from electrode with max SA to closest SOZ
    [~,I] = max(sa);
    SALoc = locs(I,:);
    SADist = min(vecnorm(SALoc - locs(soz,:),2,2));
    allSADist = [allSADist,SADist];
    
    % Distance from electrode with max seq freq of interest to its closest SOZ
    seq_freq = sum(~isnan(seq_matrix),2);
    [~,most_freq] = max(seq_freq);
    freqLoc = locs(most_freq,:);
    freqDist = min(vecnorm(freqLoc - locs(soz,:),2,2));
    allFreqDist = [allFreqDist,freqDist];
    
    % Distance from every electrode to its closest SOZ
    allLocs = zeros(size(locs,1),1);
    for i = 1:length(allLocs)
        allLocs(i) = min(vecnorm(locs(i,:) - locs(soz,:),2,2));
    end
    allAllDist =[allAllDist;allLocs];
    
    % Distance from electrodes with spikes to closest SOZ electrode
    spikeLocs = locs(ch_w_spikes,:);
    spikeDist = zeros(size(spikeLocs,1),1);
    for i = 1:size(spikeLocs,1)
        spikeDist(i) = min(vecnorm(spikeLocs(i,:)-locs(soz,:),2,2));
    end
    allSpikeDist = [allSpikeDist;spikeDist];
    
    % Get distance between electrode with max SA and electrode with max seq
    % freq
    
    allSAToFreqDist = [allSAToFreqDist;vecnorm(freqLoc-SALoc)];
    
    if whichPt == 8, exampleSA = sa; end
    if whichPt == 8, exampleChInfluence = chInfluence; end
    
    %{
    
    if doPlots == 2
        % Pretty plot
        fig = figure;
        circSize = 300;
        set(gcf,'Position',[200 100 1200 600])
        [ha, pos] = tight_subplot(1,2,[0 .03],[.05 .08],[.05 .01]); 
        % Plot connections and surface area for biggest SA channel
        axes(ha(1));
        scatter3(locs(:,1),locs(:,2),locs(:,3),circSize,'k','linewidth',2);
        hold on
        [~,I] = max(sa);
        scatter3(locs(I,1),locs(I,2),locs(I,3),circSize,'g','filled');
        downstream = chInfluence{I};
        
        for j = 1:length(downstream)
            dp = locs(downstream(j),:) - locs(I,:);
            quiver3(locs(I,1),locs(I,2),locs(I,3),dp(1), dp(2), dp(3),...
                'color','k','linewidth',2,'maxheadsize',0.4);
        end
        scatter3(locs(I,1),locs(I,2),locs(I,3),circSize,'g','filled');

        downstream = [I,downstream];
        down_locs = locs(downstream,:);
        
        tri = delaunay(down_locs(:,1),down_locs(:,2));
        %{
        cv = trisurf(tri,down_locs(:,1),down_locs(:,2),...
            down_locs(:,3),'facecolor','r');
        alpha(cv,0.05) 
        %}
        title('Downstream electrode connections')
        xticklabels([])
        yticklabels([])
        zticklabels([])
        set(gca,'fontsize',25)
        view([-0.5 -0.5 0.2])
        annotation('textbox',[0.03 0.73 0.2 0.2],'String','A','EdgeColor','none','fontsize',30);
        
        
        axes(ha(2));
        cv = trisurf(tri,down_locs(:,1),down_locs(:,2),...
            down_locs(:,3),'facecolor','r');
        alpha(cv,0.2) 
        hold on
        scatter3(locs(:,1),locs(:,2),locs(:,3),circSize,'k','linewidth',2);
        hold on
        [~,I] = max(sa);
        scatter3(locs(I,1),locs(I,2),locs(I,3),circSize,'g','filled');
        
        title('Area connecting downstream electrodes');
        set(gca,'fontsize',25)
        xticklabels([])
        yticklabels([])
        zticklabels([])
        view([-0.5 -0.5 0.2])
        annotation('textbox',[0.52 0.73 0.2 0.2],'String','B','EdgeColor','none','fontsize',30);
        fig.GraphicsSmoothing = 'off'; 
        
        %f2 = myaa(2);
        %pause
        print(fig,[saveFolder,'influence_pretty_',sprintf('%s',pt(whichPt).name)],'-dpng');
        %close(f2)
        close(fig)
    
        
    
    %% Make plots
    elseif doPlots == 1
        figure
        set(gcf,'Position',[200 100 1400 300])
        
        % Plot connections and surface area for biggest SA channel
        subplot(1,3,1)
        scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
        hold on
        [~,I] = max(sa);
        scatter3(locs(I,1),locs(I,2),locs(I,3),100,'g','filled');
        downstream = chInfluence{I};
        %{
        for j = 1:length(downstream)
            dp = locs(downstream(j),:) - locs(I,:);
            quiver3(locs(I,1),locs(I,2),locs(I,3),dp(1), dp(2), dp(3),...
                'color','k');
        end
        %}
        downstream = [I,downstream];
        down_locs = locs(downstream,:);
        
        tri = delaunay(down_locs(:,1),down_locs(:,2));
        cv = trisurf(tri,down_locs(:,1),down_locs(:,2),...
            down_locs(:,3),'facecolor','r');
        alpha(cv,0.05) 
        title('SA of influence of channel with max area')
        
        % Plot SA of influence for all channels
        subplot(1,3,2)
        gs = fh_map(50);
        [Y,E] = discretize(sa,size(gs,1));
        scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
        hold on
        scatter3(locs(isnan(Y)==0,1),locs(isnan(Y)==0,2),locs(isnan(Y)==0,3),100,gs(Y(isnan(Y)==0),:));
        [~,I] = max(sa);
        scatter3(locs(I,1),locs(I,2),locs(I,3),100,gs(Y(I),:),'filled');
        scatter3(locs(soz,1),locs(soz,2),locs(soz,3),30,'k','filled');
        title('SA of influence')
        
        % Plot seq frequency
        gs = fh_map(50);
        [Y,E] = discretize((seq_freq),size(gs,1));
        [~,most_freq] = max(seq_freq);
        
        subplot(1,3,3)
        scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
        hold on
        scatter3(locs(seq_freq~=0,1),locs(seq_freq~=0,2),locs(seq_freq~=0,3),...
            100,gs(Y(seq_freq~=0),:));
        scatter3(locs(most_freq,1),locs(most_freq,2),locs(most_freq,3),...
            100,gs(Y(most_freq),:),'filled');
        scatter3(locs(soz,1),locs(soz,2),locs(soz,3),30,'k','filled');
        title('Sequence frequency');
        
        pause
        print(gcf,[saveFolder,'influence_',sprintf('%s',pt(whichPt).name)],'-dpng');
        close(gcf)
    

    end
    %}
    
end

%{
figure
subplot(1,3,1)
scatter(allSF,allAllDist);
subplot(1,3,2)
scatter(allDegPref(allSpikers==1),allAllDist(allSpikers==1));
subplot(1,3,3)
scatter(allSAs(allSpikers==1),allAllDist(allSpikers==1));
[rho,pval] = corr(allSAs(allSpikers==1),allAllDist(allSpikers==1),'Type','Spearman')
%}

%% Does SA do better than chance?
[pFreqSA,h3,stats3] = ranksum(allFreqDist,allSADist);
%pFreqSA
[pAllSA,h4,stats4] = ranksum(allSADist,allAllDist);
%pAllSA
%[pSpikeSA,h5,stats5] = ranksum(allSADist,allSpikeDist);
[pFreqAll,~,~] = ranksum(allFreqDist,allAllDist);



%% Do Plots

%% Parameters and plot initialization
fontsizes = 15;

figure
set(gcf,'Position',[270 12 817 793]);
[ha,pos] = tight_subplot(3,2,[0.09 0.03],[0.03 0.03],[0.07 0.02]);
set(ha(5),'Position',[pos{5}(1) pos{5}(2) pos{5}(3)*2 pos{5}(4)]);
delete(ha(6));

%% Plot downstream connections
axes(ha(1))
imagesc(exampleChCh)
colorbar
title('Number of downstream spike connections');
xlabel('Downstream electrode #');
ylabel('Leading electrode #');
set(gca,'FontSize',fontsizes)
annotation('textbox',[0.02 0.8 0.2 0.2],'String','A','EdgeColor','none','fontsize',25);

%% Plot significant downstream connections
axes(ha(2))
imagesc(exampleChCh > exampleX)
colormap(ha(2),flipud(gray));
title('Frequent downstream spike connections');
xlabel('Downstream electrode #');
%ylabel('Leading electrode #');
set(gca,'FontSize',fontsizes)
yticklabels([])
annotation('textbox',[0.51 0.8 0.2 0.2],'String','B','EdgeColor','none','fontsize',25);

%% Plot downstream connections for single electrode
axes(ha(3));
circSize = 150;
locs = pt(8).electrodeData.locs(:,2:4);
scatter3(locs(:,1),locs(:,2),locs(:,3),circSize,'k','linewidth',2);
hold on
[~,I] = max(exampleSA);
scatter3(locs(I,1),locs(I,2),locs(I,3),circSize,'g','filled');
downstream = exampleChInfluence{I};

for j = 1:length(downstream)
    dp = locs(downstream(j),:) - locs(I,:);
    quiver3(locs(I,1),locs(I,2),locs(I,3),dp(1), dp(2), dp(3),...
        'color','k','linewidth',2,'maxheadsize',0.4);
end
scatter3(locs(I,1),locs(I,2),locs(I,3),circSize,'g','filled');
xticklabels([])
yticklabels([])
zticklabels([])
view(118.1000,2.8000);
title(sprintf('Downstream electrodes for electrode %d',I));
set(gca,'FontSize',fontsizes)
annotation('textbox',[0.02 0.46 0.2 0.2],'String','C','EdgeColor','none','fontsize',25);

%% Plot area of influence for a single electrode
axes(ha(4));
downstream = [I,downstream];
down_locs = locs(downstream,:);
tri = delaunay(down_locs(:,1),down_locs(:,2));
cv = trisurf(tri,down_locs(:,1),down_locs(:,2),...
    down_locs(:,3),'facecolor','r');
alpha(cv,0.2) 
hold on
scatter3(locs(:,1),locs(:,2),locs(:,3),circSize,'k','linewidth',2);
hold on
[~,I] = max(exampleSA);
scatter3(locs(I,1),locs(I,2),locs(I,3),circSize,'g','filled');
xticklabels([])
yticklabels([])
zticklabels([])
title(sprintf('Area of influence for electrode %d',I));
view(118.1000,2.8000);
set(gca,'FontSize',fontsizes)
annotation('textbox',[0.51 0.46 0.2 0.2],'String','D','EdgeColor','none','fontsize',25);

%% Plot bar graph showing overall performance
axes(ha(5));
prices = [mean(allAllDist) mean(allSADist) mean(allFreqDist)];
bar(prices)
title(sprintf(['Average distance across patients from electrode of interest\n to closest ',...
    'seizure onset zone electrode']))
ylabel(sprintf('Average distance (mm)'));

xticklabels({'All electrodes','Max area of influence','Max spike frequency'})
set(gca,'FontSize',fontsizes)
%fix_xticklabels(gca,0.1,{'FontSize',fontsizes});

% Plot p-values
if pFreqSA < 0.001
    textFreqSA = 'p < 0.001';
else
    textFreqSA = sprintf('p = %1.3f',pFreqSA);
end
hold on
plot([2.1 3], [max(prices)+1 max(prices)+1],'k')
text(2.5,max(prices)+4,textFreqSA,'HorizontalAlignment','center',...
        'fontsize',fontsizes);
    
if pAllSA < 0.001
    textAllSA = 'p < 0.001';
else
    textAllSA = sprintf('p = %1.3f',pAllSA);
end


hold on
plot([1 1.9], [max(prices)+1 max(prices)+1],'k')
text(1.5,max(prices)+4,textAllSA,'HorizontalAlignment','center',...
        'fontsize',fontsizes);
    
if pFreqAll < 0.001
    textFreqAll = 'p < 0.001';
else
    textFreqAll = sprintf('p = %1.3f',pFreqAll);
end
plot([1 3], [max(prices)+8 max(prices)+8],'k');
text(2,max(prices)+11,textFreqAll,'HorizontalAlignment','center',...
        'fontsize',fontsizes);
    
ylim([0 max(prices) + 15]);
annotation('textbox',[0.02 0.11 0.2 0.2],'String','E','EdgeColor','none','fontsize',25);
%pause
print(gcf,[destFolder,'influence_nonAA'],'-depsc');
eps2pdf([destFolder,'influence_nonAA','.eps'])


f2= myaa(2);
%pause
print(f2,[destFolder,'influence_AA'],'-dpng');
%eps2pdf([destFolder,'influence_AA','.eps'])

fprintf('The average largest area of influence was %1.1f cm squared (range %1.1f-%1.1f).\n',...
    mean(allMaxSAs/100),min(allMaxSAs/100),max(allMaxSAs/100));

fprintf('The average distance between biggest SA electrode and nearest SOZ electrode is:\n%1.1f (range %1.1f-%1.1f)\n',...
    mean(allSADist), min(allSADist), max(allSADist));

fprintf(['The average distance between the electrodes with max SA and\n'...
    'the electrodes with max frequency is: %1.1f (range %1.1f-%1.1f)\n and there were %d of %d ',...
    'patients where these electrodes were the same\n'],mean(allSAToFreqDist),...
    min(allSAToFreqDist),max(allSAToFreqDist),sum((allSAToFreqDist==0)),length(allSAToFreqDist==0));



end