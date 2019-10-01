function CInfluence(pt,cluster,whichPts)

%{

My area of influence analysis.

this take a pt struct, a cluster struct, and calculates
distances of various electrodes of interest, including the area with the
largest area of influence, from the nearest SOZ.

%}

% Parameters
doPoster = 1;
doPlots = 0; %0 = no, 1=normal, 2=pretty
plotConn = 0;
removeTies = 1;
doBootstrap = 0; % 0 = no, 1 = standard, 2 = permute rows, 3 = permute columns
alpha1 = 95;
map_text = 'jet';
fh_map = str2func(map_text);
new_bootstrap = 0;

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
outcome_all = [];
temp_lobe_all = [];
allMeanDist = [];
sa_resected = [];
sf_resected = [];
sa_sf_corr = [];

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
   % error('Warning, not doing correct patients!\n');
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
    
    % outcomes
    %outcome = getOutcome(pt,whichPt);
    outcome = str2num(pt(whichPt).clinical.outcome(end));
    outcome_all = [outcome_all;outcome];
    
    % SOZ
    szOnsetText = pt(whichPt).clinical.seizureOnset;
    if contains(szOnsetText,'TL') == 1
        tempLobe = 1;
    else
        tempLobe = 0;
    end
    temp_lobe_all = [temp_lobe_all;tempLobe];
    
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
    
    spike_count_per_ch = sum(chCh,1);
    
    if whichPt == 8, exampleChCh = chCh; end
    
    
    
    %% Get significant connections
    if doBootstrap == 2
        % This is a new analysis where I compare the number of connections
        % to that of a matrix where I randomly permute just the rows of the
        % connection matrix. In this random matrix, the number of times
        % downstream electrodes are activated is held constant, but this
        % set of downstream connections is attributed to a random leader
        % electrode.
        
        nboot = 1e3;
        boot_con = zeros(nboot,nchs,nchs);
        for ib = 1:nboot
            
            % Get random permutation from 1:nchs
            p = randperm(nchs);
            
            % shuffle the connection matrix according to this permutation
            temp_chCh =  chCh(p,:);
            boot_con(ib,:,:) = temp_chCh;
     
        end
        
        % Now go through each row, and then for each column in that row,
        % find the 95% number of connections in the bootstrap matrix.
        % This is the minimum connection number needed to achieve
        % significance.
        
        sig_con = zeros(nchs,nchs);
        chInfluence = cell(nchs,1);
        for i = 1:nchs
            for j = 1:nchs
                
                % sort the bootstrap connections for that element
                boot_element = sort(boot_con(:,i,j));
                perc = prctile(boot_element,alpha1);
                
                % See if the number of connections is higher than this
                if chCh(i,j) > perc
                    sig_con(i,j) = 1;
                    chInfluence{i} = [chInfluence{i},j];
                end
                
            end
        end
        
    elseif doBootstrap == 3
        
        % Now permute the columns (so each leader electrode will keep the
        % same general downstream number of connections, but they will be
        % assigned to a random permutation of downstream connections).
        
        nboot = 1e3;
        boot_con = zeros(nboot,nchs,nchs);
        for ib = 1:nboot
            
            % Get random permutation from 1:nchs
            p = randperm(nchs);
            
            % shuffle the connection matrix according to this permutation
            temp_chCh =  chCh(:,p);
            boot_con(ib,:,:) = temp_chCh;
     
        end
        
        % Now go through each row, and then for each column in that row,
        % find the 95% number of connections in the bootstrap matrix.
        % This is the minimum connection number needed to achieve
        % significance.
        
        sig_con = zeros(nchs,nchs);
        chInfluence = cell(nchs,1);
        for i = 1:nchs
            for j = 1:nchs
                
                % sort the bootstrap connections for that element
                boot_element = sort(boot_con(:,i,j));
                perc = prctile(boot_element,alpha1);
                
                % See if the number of connections is higher than this
                if chCh(i,j) > perc
                    sig_con(i,j) = 1;
                    chInfluence{i} = [chInfluence{i},j];
                end
                
            end
        end
    
    elseif doBootstrap == 1
        
        % Here, for each permutation, I am constructing a chCh matrix where
        % I am distributing the true total number of connections randomly
        % across all elements of the nch by nch matrix.
        
        ncons = sum(sum(chCh));
        nboot = 1e3;
        max_size = nchs*nchs;
        chCh_all = zeros(nboot,nchs,nchs);
        chCh_diff_all = zeros(nboot,nchs,nchs);
        chCh_diff = zeros(nchs,nchs);
        for ib = 1:nboot
            ib
            if mod(ib,100) == 0
                fprintf('Doing %d of %d\n', ib,nboot);
            end
            chCh_f = zeros(nchs,nchs);
            if new_bootstrap == 0
                for j = 1:ncons
                    chCh_f(randi(max_size)) = chCh_f(randi(max_size)) + 1;

                end
                chCh_all(ib,:,:) = chCh_f;
            elseif new_bootstrap == 1
                % weight the probability of getting a connection by the
                % spike rate
               
                    
                rows = randsample(nchs,ncons,true,spike_count_per_ch);
                cols = zeros(ncons,1);
                for j = 1:ncons
                    not_allowed = rows(j); % spike can't travel from one channel to same channel
                    new_weights = spike_count_per_ch; new_weights(not_allowed) = 0;
                    cols(j) = randsample(nchs,1,true,new_weights);  
                end
                
                for j = 1:ncons
                   chCh_all(ib,rows(j),cols(j)) = chCh_all(ib,rows(j),cols(j)) + 1;
                end
                
               
            end
            
            % Also calculate the difference between i,j and j,i
            for i = 1:size(chCh_f,1)
                for j = 1:size(chCh_f,2)
                    chCh_diff(i,j) = chCh_f(i,j) - chCh_f(j,i);
                end
            end
            chCh_diff_all(ib,:,:) = chCh_diff;
            
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
    
        % Get the 95% of differences 
        s_diff = sort(chCh_diff_all(:));
        perc_diff = prctile(s_diff,alpha1);
        
    end
    
    n_spikes_ch = sum(~isnan(seq_matrix),2);
    
    if doBootstrap < 2
    % Assume poisson distribution (produces same result as permutation
    % test)
    ncons = sum(sum(chCh));
    lambda = ncons/nchs^2;
    X = poissinv(alpha1/100,lambda);
    minCount = X;
    fprintf(['By poisson assumption, the number of counts is:\n'...
        '%d\n\n'],X);
    
    
    %if perc~=X, error('What\n'); end
    
    if whichPt == 8, exampleX = X; end
    
    %% Now find connections that are more frequent by chance
    n_spikes = length(all_spike_times);
    lambda_spikes = n_spikes/nchs;
    X_spikes = poissinv(alpha1/100,lambda_spikes);
    minCountSpikes = X_spikes;
    
    
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
        %sa(i) = areaConnecting(down_locs,0);
    end
    
    
    %% How far is the channel with max area of influence from SOZ?
    allSF = [allSF;n_spikes_ch];
    allSAs = [allSAs;sa];
    allMaxSAs = [allMaxSAs;(max(sa))];
    
    %% Make AAN plot
    if 0
    figure
    scatter3(locs(:,1),locs(:,2),locs(:,3),200,'k');
    %parula_c = parula(2);
    parula_c = [0 0.4470 0.7410;
    0.8500 0.3250 0.0980];
    hold on
    seq_freq = sum(~isnan(seq_matrix),2);
    [~,temp_max_sf] = max(seq_freq);
    
    h_soz = scatter3(locs(soz,1),locs(soz,2),locs(soz,3),200,parula_c(1,:),'filled');
    h_sf = scatter3(locs(temp_max_sf,1),locs(temp_max_sf,2),...
        locs(temp_max_sf,3),100,parula_c(2,:),'filled');
    legend([h_soz,h_sf],{'Seizure onset zone','Max spike frequency'},'location',...
        'northwest');
    view(-14.7,30)
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    zticklabels([])
    yticklabels([])
    xticklabels([])
    grid off
    set(gca,'fontsize',25)
    print(gcf,[destFolder,'brain_ex_AAN'],'-depsc');
    end
    

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
    allAllDist =[allAllDist;mean(allLocs)];
    
    % Distance from electrodes with spikes to closest SOZ electrode
    if 0 
    spikeLocs = locs(ch_w_spikes,:);
    spikeDist = zeros(size(spikeLocs,1),1);
    for i = 1:size(spikeLocs,1)
        spikeDist(i) = min(vecnorm(spikeLocs(i,:)-locs(soz,:),2,2));
    end
    allSpikeDist = [allSpikeDist;spikeDist];
    end
    
    % Get distance between electrode with max SA and electrode with max seq
    % freq
    
    allSAToFreqDist = [allSAToFreqDist;vecnorm(freqLoc-SALoc)];
    allMeanDist = [allMeanDist;mean(allLocs)];
    
    if whichPt == 8, exampleSA = sa; end
    if whichPt == 8, exampleSF = seq_freq; end
    if whichPt == 8, exampleChInfluence = chInfluence; end
    
     % Get list of resected electrodes
    resec_elecs = pt(whichPt).resecElecs;
    
    %% Was max SA electrode resected?
    % Get electrode with max sa
    [~,sa_max_ident] = max(sa);
    sa_resected = [sa_resected;ismember(sa_max_ident,resec_elecs)];
    
    %% Was max SF electrode resected?
    % Get electrode with max SF
    [~,sf_max_ident] = max(seq_freq);
    sf_resected = [sf_resected;ismember(sf_max_ident,resec_elecs)];
    
    %% Calculate spearman rank correlation between SF and SA
    sa_sf_corr = [sa_sf_corr;corr(seq_freq,sa,'Type','Spearman')];
    
    %% Example plot to convince myself I am doing the resection analysis correctly
%{
    figure
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k','linewidth',2);
    hold on
    scatter3(locs(sa_max_ident,1),locs(sa_max_ident,2),locs(sa_max_ident,3),100,'b','filled');
    scatter3(locs(sf_max_ident,1),locs(sf_max_ident,2),locs(sf_max_ident,3),100,'r','filled');
    scatter3(locs(resec_elecs,1),locs(resec_elecs,2),locs(resec_elecs,3),20,'k','filled');
    fprintf('Was max SA resected? %d\n',ismember(sa_max_ident,resec_elecs));
    fprintf('Was max SF resected? %d\n\n\n',ismember(sf_max_ident,resec_elecs));
    pause
    close(gcf)
    
    %}
    
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

aes_plot_aoi(allAllDist,allSADist,allFreqDist)

%% Make AAN plots
figure
np = size(allAllDist,1);
%{
scatter(ones(np,1) + ones(np,1).*randn(np,1)/15,allFreqDist,100,'filled');
hold on
scatter(2*ones(np,1) + ones(np,1).*randn(np,1)/15,allAllDist,100,'filled');
%}
h = boxplot([allFreqDist',allAllDist]);
h2 = h(:);
for ih = 1:length(h2)
    set(h2(ih,:),'LineWidth',2);
end
hold on
[pFreqAll,h5,stats5] = signrank(allFreqDist,allAllDist);
if pFreqAll < 0.001
    textFreqAll = 'p < 0.001***';
elseif pFreqAll < 0.01
    textFreqAll = sprintf('p = %1.3f** (Wilcoxon signed-rank)',pFreqAll);
elseif pFreqAll < 0.05
    textFreqAll = sprintf('p = %1.3f*',pFreqAll);
else
    textFreqAll = sprintf('p = %1.3f',pFreqAll);
end
max_point = max([allFreqDist,allAllDist']);
plot([1 2],[max_point+2 max_point+2],'k','linewidth',2)
text(1.5,max_point+7,textFreqAll,'HorizontalAlignment','Center','FontSize',20)
xticks([1 2])
xticklabels({'Max spike rate','All electrodes'})
ylabel({'Distance from nearest','seizure onset zone (mm)'})
set(gca,'fontsize',20)
set(gca,'xlim',[0.7 2.3])
set(gca,'ylim',[0 max_point+10])
print([destFolder,'aan_soz'],'-depsc')

%% Get average SRC between SA and SF
% I'll just do a plain average across patients since I am not doing
% significance testing
rho_sa_sf_avg = mean(sa_sf_corr);

%% Compare clinical outcome for pts with resected max SA and those without
[p_sa_resec,h_sa_resec,sa_stats_resec] = ...
    ranksum(outcome_all(sa_resected == 1),outcome_all(sa_resected == 0));

fprintf(['P value for different outcome between resected max area of influence'...
    'and not is:\n%1.1e, ranksum = %1.1f\n'],...
    p_sa_resec,sa_stats_resec.ranksum);

ilae_sa_resec = outcome_all(sa_resected == 1);%getILAE(outcome_all(sa_resected == 1));
ilae_sa_noresec = outcome_all(sa_resected == 0);%getILAE(outcome_all(sa_resected == 0));

fprintf('Median ilae for sa resec = %1.1f, for sa no resec %1.1f\n',...
    median(ilae_sa_resec), median(ilae_sa_noresec));


[~,~, u_mat] = ranksum_erin(outcome_all(sa_resected == 1),outcome_all(sa_resected == 0));
[~,~, u_erin] = ranksum_erin(outcome_all(sa_resected == 1),outcome_all(sa_resected == 0));
fprintf('Ranksum matlab  = %1.1f, ranksum erin = %1.1f\n',u_mat,u_erin);

%% Compare clinical outcome for pts with resected max SF and those without
[p_sf_resec,h_sf_resec,stats_sf_resec] = ...
    ranksum(outcome_all(sf_resected == 1),outcome_all(sf_resected == 0));

fprintf(['P value for different outcome between resected max spike frequency'...
    'and not is:\n%1.1e, ranksum = %1.1f\n'],...
    p_sf_resec,stats_sf_resec.ranksum);

ilae_sf_resec = outcome_all(sf_resected == 1);%getILAE(outcome_all(sf_resected == 1));
ilae_sf_noresec = outcome_all(sf_resected == 0);%getILAE(outcome_all(sf_resected == 0));

fprintf('Median ilae for sf resec = %1.1f, for sf no resec %1.1f\n',...
    median(ilae_sf_resec), median(ilae_sf_noresec));

[~,~, u_mat] = ranksum_erin(outcome_all(sf_resected == 1),outcome_all(sf_resected == 0));
[u_erin] = getStandardStats(outcome_all(sf_resected == 1),outcome_all(sf_resected == 0),'rs');
fprintf('Ranksum matlab  = %1.1f, ranksum erin = %1.1f\n',u_mat,u_erin);

%% Does SA do better than chance?
%{
[pFreqSA,h3,stats3] = ranksum(allFreqDist,allSADist);
%pFreqSA
[pAllSA,h4,stats4] = ranksum(allSADist,allAllDist);
%pAllSA
%[pSpikeSA,h5,stats5] = ranksum(allSADist,allSpikeDist);
[pFreqAll,h5,stats5] = ranksum(allFreqDist,allAllDist);
%}

%% SA vs freq
[pFreqSA,h3,stats3] = signrank(allFreqDist,allSADist);

fprintf('P-value for max freq vs max SA is %1.1e, signed-rank = %1.1f\n',...
    pFreqSA,stats3.signedrank);

w_mat = signrank_erin(allFreqDist',allSADist');
w_erin = getStandardStats(allFreqDist',allSADist','sr');
fprintf('Signrank matlab  = %1.1f, signrank erin = %1.1f\n',w_mat,w_erin);


%% SA vs all
[pAllSA,h4,stats4] = signrank(allSADist',allAllDist);

fprintf('P-value for max SA vs all is %1.1e, signed-rank = %1.1f\n',...
    pAllSA,stats4.signedrank);

w_mat = signrank_erin(allSADist',allAllDist);
w_erin = getStandardStats(allSADist',allAllDist,'sr');
fprintf('Signrank matlab  = %1.1f, signrank erin = %1.1f\n',w_mat,w_erin);

%% Freq vs all

[pFreqAll,h5,stats5] = signrank(allFreqDist,allAllDist);
fprintf('P-value for max freq vs all is %1.1e, signed-rank = %1.1f\n',...
    pFreqAll,stats5.signedrank);
w_mat = signrank_erin(allFreqDist',allAllDist);
w_erin = getStandardStats(allFreqDist',allAllDist,'sr');
fprintf('Signrank matlab  = %1.1f, signrank erin = %1.1f\n',w_mat,w_erin);



%% Do Plots

%{
%% First, plot for poster
fontsizes = 19;
figure
set(gcf,'Position',[71 107 1370 695]);
[ha,pos] = tight_subplot(2,3,[0.1 0.03],[0.05 0.05],[0.07 0.02]);
delete(ha(4));
delete(ha(5));
set(ha(6),'Position',[pos{6}(1)-0.07 pos{6}(2) pos{6}(3)+0.07 pos{6}(4)]);

axes(ha(1))
imagesc(exampleChCh)
colorbar
title('Downstream spike connections');
xlabel('Downstream electrode #');
ylabel('Leading electrode #');
set(gca,'FontSize',fontsizes)

axes(ha(2))
imagesc(exampleChCh > exampleX)
colormap(ha(2),flipud(gray));
title('Frequent connections');
xlabel('Downstream electrode #');
%ylabel('Leading electrode #');
set(gca,'FontSize',fontsizes)
yticklabels([])

axes(ha(3));
circSize = 150;
locs = pt(8).electrodeData.locs(:,2:4);
[~,I] = max(exampleSA);
downstream = exampleChInfluence{I};
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

axes(ha(6));
prices = [mean(allAllDist) mean(allSADist) mean(allFreqDist)];
bar(prices)
title('Distance to seizure onset zone')
ylabel(sprintf('Distance to closest\nseizure onset zone electrode\(mm)'));

xticklabels({'All electrodes','Max area of influence','Max spike rate'})
set(gca,'FontSize',fontsizes)
set(gca,'xlim',[0.5 3.5]);
%fix_xticklabels(gca,0.1,{'FontSize',fontsizes});

% Plot p-values
if pFreqSA < 0.001
    textFreqSA = 'p < 0.001***';
else
    textFreqSA = sprintf('p = %1.3f',pFreqSA);
end
hold on
plot([2.1 3], [max(prices)+1 max(prices)+1],'k')
text(2.5,max(prices)+4,textFreqSA,'HorizontalAlignment','center',...
        'fontsize',fontsizes);
    
if pAllSA < 0.001
    textAllSA = 'p < 0.001***';
elseif pAllSA < 0.01
    textAllSA = sprintf('p = %1.3f**',pAllSA);
elseif pAllSA < 0.05
    textAllSA = sprintf('p = %1.3f*',pAllSA);
else
    textAllSA = sprintf('p = %1.3f',pAllSA);
end


hold on
plot([1 1.9], [max(prices)+1 max(prices)+1],'k')
text(1.5,max(prices)+4,textAllSA,'HorizontalAlignment','center',...
        'fontsize',fontsizes);
    
if pFreqAll < 0.001
    textFreqAll = 'p < 0.001***';
elseif pFreqAll < 0.01
    textFreqAll = sprintf('p = %1.3f**',pFreqAll);
elseif pFreqAll < 0.05
    textFreqAll = sprintf('p = %1.3f*',pFreqAll);
else
    textFreqAll = sprintf('p = %1.3f',pFreqAll);
end
plot([1 3], [max(prices)+8 max(prices)+8],'k');
text(2,max(prices)+11,textFreqAll,'HorizontalAlignment','center',...
        'fontsize',fontsizes);
    
ylim([0 max(prices) + 15]);
%pause
print(gcf,[destFolder,'influence_nonAA_poster'],'-depsc');
eps2pdf([destFolder,'influence_nonAA_poster','.eps'])
%}


%% Parameters and plot initialization
fontsizes = 20;

figure
if doPoster == 1
    set(gcf,'Position',[107 12 1200 793]);
    [ha,pos] = tight_subplot(3,3,[0.11 0.05],[0.04 0.04],[0.06 0.04]);
    set(ha(7),'Position',[pos{7}(1) pos{7}(2) pos{7}(3)*1.5 pos{7}(4)]);
else
    set(gcf,'Position',[107 12 946 793]);
    [ha,pos] = tight_subplot(3,3,[0.09 0.03],[0.03 0.03],[0.07 0.02]);
    set(ha(7),'Position',[pos{7}(1) pos{7}(2) pos{7}(3)*3 pos{7}(4)]);
end


delete(ha(8));
delete(ha(9));

%% Plot downstream connections
axes(ha(1))
imagesc(exampleChCh)
colorbar
title('HUP078 spike connections');
xlabel('Downstream electrode');
ylabel('Leading electrode');
set(gca,'FontSize',fontsizes)
if doPoster == 0
    annotation('textbox',[0.02 0.785 0.2 0.2],'String','A','EdgeColor','none','fontsize',25);
end

%% Plot significant downstream connections
axes(ha(2))
imagesc(exampleChCh > exampleX)
colormap(ha(2),flipud(gray));
title('Frequent connections');
xlabel('Downstream electrode');
%ylabel('Leading electrode #');
set(gca,'FontSize',fontsizes)
yticklabels([])
if doPoster == 0
annotation('textbox',[0.35 0.785 0.2 0.2],'String','B','EdgeColor','none','fontsize',25);
end

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
title(sprintf('Electrode %d''s downstream electrodes',I));
set(gca,'FontSize',fontsizes)
if doPoster == 0
    annotation('textbox',[0.67 0.785 0.2 0.2],'String','C','EdgeColor','none','fontsize',25);
end

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
if doPoster == 0
    annotation('textbox',[0.02 0.45 0.2 0.2],'String','D','EdgeColor','none','fontsize',25);
end

%% Plot area of influence of all electrodes
axes(ha(5));
scatter3(locs(:,1),locs(:,2),locs(:,3),circSize,'k','linewidth',2);
hold on
scatter3(locs(:,1),locs(:,2),locs(:,3),circSize,exampleSA,'filled')
xticklabels([])
yticklabels([])
zticklabels([])
title(sprintf('Area of influence, all electrodes'));
view(118.1000,2.8000);
set(gca,'FontSize',fontsizes)
if doPoster == 0
    annotation('textbox',[0.35 0.45 0.2 0.2],'String','E','EdgeColor','none','fontsize',25);
end

%% Plot spike frequency of all electrodes
axes(ha(6));
scatter3(locs(:,1),locs(:,2),locs(:,3),circSize,'k','linewidth',2);
hold on
scatter3(locs(:,1),locs(:,2),locs(:,3),circSize,exampleSF,'filled')
xticklabels([])
yticklabels([])
zticklabels([])
title(sprintf('Spike frequency, all electrodes'));
view(118.1000,2.8000);
set(gca,'FontSize',fontsizes)
if doPoster == 0
    annotation('textbox',[0.67 0.45 0.2 0.2],'String','F','EdgeColor','none','fontsize',25);
end


%% Plot bar graph showing overall performance
axes(ha(7));
prices = [mean(allAllDist) mean(allSADist) mean(allFreqDist)];
bar(prices)
title(sprintf(['All patients: Average distance to closest\n',...
    'seizure onset zone electrode']))
ylabel(sprintf('Average distance (mm)'));
set(gca,'xlim',[0.7 3.3]);
xticklabels({'All electrodes','Max area of influence','Max spike rate'})
set(gca,'FontSize',fontsizes)
%fix_xticklabels(gca,0.1,{'FontSize',fontsizes});

% Plot p-values
if pFreqSA < 0.001
    textFreqSA = 'p < 0.001***';
else
    textFreqSA = sprintf('p = %1.3f',pFreqSA);
end
hold on
plot([2.1 3], [max(prices)+1 max(prices)+1],'k')
text(2.5,max(prices)+4,textFreqSA,'HorizontalAlignment','center',...
        'fontsize',fontsizes);
    
if pAllSA < 0.001
    textAllSA = 'p < 0.001***';
elseif pAllSA < 0.01
    textAllSA = sprintf('p = %1.3f**',pAllSA);
elseif pAllSA < 0.05
    textAllSA = sprintf('p = %1.3f*',pAllSA);
else
    textAllSA = sprintf('p = %1.3f',pAllSA);
end


hold on
plot([1 1.9], [max(prices)+1 max(prices)+1],'k')
text(1.5,max(prices)+4,textAllSA,'HorizontalAlignment','center',...
        'fontsize',fontsizes);
    
if pFreqAll < 0.001
    textFreqAll = 'p < 0.001***';
elseif pFreqAll < 0.01
    textFreqAll = sprintf('p = %1.3f**',pFreqAll);
elseif pFreqAll < 0.05
    textFreqAll = sprintf('p = %1.3f*',pFreqAll);
else
    textFreqAll = sprintf('p = %1.3f',pFreqAll);
end
plot([1 3], [max(prices)+8 max(prices)+8],'k');
text(2,max(prices)+11,textFreqAll,'HorizontalAlignment','center',...
        'fontsize',fontsizes);

ylim([0 max(prices) + 15]);
if doPoster == 0
annotation('textbox',[0.02 0.11 0.2 0.2],'String','G','EdgeColor','none','fontsize',25);
end
%pause
if doPoster == 1
    print(gcf,[destFolder,'influence_nonAA_poster'],'-depsc');
else
    print(gcf,[destFolder,'influence_nonAA'],'-depsc');
end
%eps2pdf([destFolder,'influence_nonAA','.eps'])


f2= myaa(2);
%pause
%print(f2,[destFolder,'influence_AA'],'-dpng');
%eps2pdf([destFolder,'influence_AA','.eps'])

fprintf('The average largest area of influence was %1.1f cm squared (range %1.1f-%1.1f).\n',...
    mean(allMaxSAs/100),min(allMaxSAs/100),max(allMaxSAs/100));

fprintf('The average distance between biggest SA electrode and nearest SOZ electrode is:\n%1.1f (range %1.1f-%1.1f)\n',...
    mean(allSADist), min(allSADist), max(allSADist));

fprintf('The average distance between highest SF electrode and nearest SOZ electrode is:\n%1.1f\n',...
    mean(allFreqDist));

fprintf('The average distance between every electrode and nearest SOZ electrode is:\n%1.1f\n',...
    mean(allAllDist));

fprintf(['The average distance between the electrodes with max SA and\n'...
    'the electrodes with max frequency is: %1.1f (range %1.1f-%1.1f)\n and there were %d of %d ',...
    'patients where these electrodes were the same\n'],mean(allSAToFreqDist),...
    min(allSAToFreqDist),max(allSAToFreqDist),sum((allSAToFreqDist==0)),length(allSAToFreqDist==0));


% Outcome correlation

good_outcome = outcome_all <= 2.25;
bad_outcome = outcome_all > 2.25;
[p,~,stats6] = ranksum((allSADist(good_outcome)),allFreqDist(good_outcome));
fprintf(['When ILAE 1-3 considered, difference between most spikey and biggest SA:\n'...
    '%1.1e\n'],p);


[p,~,stats] = ranksum(allSADist(good_outcome),allSADist(bad_outcome));
fprintf('SA dist between good outcome versus bad outcome is:\n%1.1e and ranksum %1.1f\n',...
    p,stats.ranksum);

end