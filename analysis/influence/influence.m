function influence(pt,whichPts)

% Need to decide whether to take clustering results. My intuition says no,
% but then I should have a check to throw out bad patients at least.


% Parameters
removeTies = 1;
doPlots = 0;
doBootstrap = 0;
alpha1 = 95;
map_text = 'jet';
fh_map = str2func(map_text);

[scriptFolder,resultsFolder] = fileLocations;
p1 = genpath(scriptFolder);
addpath(p1);

destFolder = [resultsFolder,'influence/'];
mkdir(destFolder)

% Initialize matrices for all patients
allFreqDist = [];
allAllDist =[];
allSADist = [];

if isempty(whichPts) == 1
    for i = 1:length(pt)
        if isempty(pt(i).seq_matrix) == 0
            whichPts = [whichPts,i];
        end
    end
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
    
    %% Get all sequences
    seq_matrix = pt(whichPt).seq_matrix;
    
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
    first_time = min(seq_matrix,[],1);
    t = (any(first_time >= (szTimes(:,1)-repmat(60,size(szTimes,1),1)) ...
        & first_time <= szTimes(:,2),2));
    seq_matrix(:,t) = [];
    fprintf('Removed %d ictal spikes \n',sum(t));
    
    %% Remove ties
    
    %% Construct a matrix of channel connections
    chCh = zeros(nchs,nchs);
    for i = 1:size(seq_matrix,2)
        seq = seq_matrix(:,i);
        spike_chs = find(isnan(seq) == 0);
        spike_times = seq((isnan(seq) == 0));
        
        % sort channels by spike time
        [~,I] = sort(spike_times);
        spike_chs = spike_chs(I);
        
        % Get all unidirectional pairwise connections that start with first
        % channel
        B = nchoosek(spike_chs,2);
        B = B(B(:,1) == spike_chs(1),:);
        
        for j = 1:size(B,1)
            chCh(B(j,1),B(j,2)) = chCh(B(j,1),B(j,2)) + 1;
        end
        
    end
    
    %% Get significant connections
    if doBootstrap == 1
        
        % Here, for each permutation, I am constructing a chCh matrix where
        % I am distributing the true total number of connections randomly
        % across all elements of the nch by nch matrix.
        ncons = sum(sum(chCh));
        nboot = 1e3;
        max_size = nchs*nchs;
        chCh_all = zeros(nboot,nchs,nchs);
        for ib = 1:nboot
            if mod(ib,100) == 0
                fprintf('Doing %d of %d\n', ib,nboot);
            end
            chCh_f = zeros(nchs,nchs);
            for j = 1:ncons
               chCh_f(randi(max_size)) = chCh_f(randi(max_size)) + 1;
            end
            chCh_all(ib,:,:) = chCh_f;
        end
        
        if 1 == 0
            figure
            imagesc(squeeze(mean(chCh_all,1)))
            colorbar
        end
        
        % Get the 95% of number of connections
        s_con = sort(chCh_all(:));
        perc = prctile(s_con,alpha1);
        
        
        fprintf(['By permutation testing, the minimum number of counts for a\n'...
        'connection to be significant is\n'...
        '%d for an alpha of %1.1fth percentile\n\n'],perc,alpha1);
        
    end
    
    % Assume poisson distribution
    ncons = sum(sum(chCh));
    lambda = ncons/nchs^2;
    X = poissinv(alpha1/100,lambda);
    minCount = X;
    fprintf(['By poisson assumption, the number of counts is:\n'...
        '%d\n\n'],X);
    
    
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
    
    %% Make plots
    if doPlots == 1
        figure
        set(gcf,'Position',[200 100 1400 300])
        
        % Plot connections and surface area for biggest SA channel
        subplot(1,3,1)
        scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
        hold on
        [~,I] = max(sa);
        scatter3(locs(I,1),locs(I,2),locs(I,3),100,'g','filled');
        downstream = chInfluence{I};
        for j = 1:length(downstream)
            dp = locs(downstream(j),:) - locs(I,:);
            quiver3(locs(I,1),locs(I,2),locs(I,3),dp(1), dp(2), dp(3),...
                'color','k');
        end
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

    end
    
end

%% Does SA do better than chance?
[pFreqSA,h3,stats3] = ranksum(allFreqDist,allSADist);
pFreqSA
[pAllSA,h4,stats4] = ranksum(allSADist,allAllDist);
pAllSA

%% Plot bars of means
figure
prices = [mean(allAllDist) mean(allSADist) mean(allFreqDist)];
bar(prices)
title(sprintf(['Average distance across patients from electrode of interest\n to closest ',...
    'seizure onset zone electrode']))
ylabel(sprintf('Average distance (mm)'));

xticklabels({'All electrodes','max SA','Electrode with max sequence frequency'})
set(gca,'FontSize',15)
fix_xticklabels(gca,0.1,{'FontSize',15});

% Plot p-values
if pFreqSA < 0.001
    textFreqSA = 'p < 0.001';
else
    textFreqSA = sprintf('p = %1.3f',pFreqSA);
end
hold on
plot([2 3], [max(prices) max(prices)],'k')
text(2.5,max(prices)+1,textFreqSA,'HorizontalAlignment','center',...
        'fontsize',15);
    
if pAllSA < 0.001
    textAllSA = 'p < 0.001';
else
    textAllSA = sprintf('p = %1.3f',pAllSA);
end
hold on
plot([1 2], [max(prices) max(prices)],'k')
text(1.5,max(prices)+1,textAllSA,'HorizontalAlignment','center',...
        'fontsize',15);



end