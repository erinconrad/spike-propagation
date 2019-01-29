function [seq_lengths,seq_times] = getSeqDist(pt,cluster,whichPts)

%{
Get the percentage of seqs that are all the same cluster

Could also get the order in which clusters get activated in sequences


%}
%% Parameters
removeTies = 1;

%% Pick Patients
if isempty(whichPts) == 1
    for i = 1:length(pt)
        if isempty(pt(i).seq_matrix) == 0
            if size(cluster(i).bad_cluster) < cluster(i).k
                whichPts = [whichPts,i];
            end
        end
    end
    
    if isequal(whichPts,[1,4,6,7,8,9,12,14,15,16,17,18,19,20,22,24,25,27,30,31]) == 0
        error('Warning, not doing correct patients!\n');
    end
end

for whichPt = whichPts
    
   %% Patient parameters
    fprintf('Doing %s\n',pt(whichPt).name);
    locs = pt(whichPt).electrodeData.locs(:,2:4);
    nchs = size(locs,1);
    szTimes = pt(whichPt).newSzTimes;
    soz = pt(whichPt).newSOZChs; 
    
    
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
    
    
    %% Get sequence lengths
    seq_lengths = zeros(size(seq_matrix,2),1);
    seq_times = zeros(size(seq_matrix,2),1);
    for s = 1:size(new_seq_matrix,2)
        currSeq = seq_matrix(:,s);
        seq_lengths(s) = sum(isnan(currSeq) == 0); 
        seq_times(s) = min(currSeq);
    end
    
    %{
   %% Get average
   run_times = pt(whichPt).runTimes;
   seq_length_bin = zeros(size(run_times,1),1);
   seq_time_bin = zeros(size(run_times,1),1);
   for i = 1:size(run_times,1)
        seq_idx = find(seq_times >= run_times(i,1) & seq_times <= run_times(i,2));
        seq_length_bin(i) = mean(seq_lengths(seq_idx));
        seq_time_bin(i) = mean(seq_times(seq_idx));
   end
   
   figure
   plot(seq_time_bin/3600,seq_length_bin);
   hold on
   for j = 1:size(szTimes,1) 
        yl = ylim;
        plot([szTimes(j,1) szTimes(j,1)]/3600,yl,'k','LineWidth',2);
   end
    
    %}
   
    
end




end