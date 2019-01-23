function allSpikeCounts(pt,cluster,whichPts)

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

all = [];
post_proc = [];
post_clust = [];
times_all_pts = [];
all_spike_num = [];

for whichPt = whichPts
    
    szTimes = pt(whichPt).newSzTimes;
    nchs = size(pt(whichPt).channels,1);
    
    %% Get spike rate before processing
    
    % Get number of spikes
    n_spikes_all = pt(whichPt).stats.nspikes;
    
    % Get all times
    all_times = sum(diff(pt(whichPt).allTimes,1,2));
    times_all_pts = [times_all_pts;all_times];
    
    % Spike rate
    spike_rate = n_spikes_all/all_times/nchs*60;
    all = [all;spike_rate];
    
    %% Get spike rate post-sequence detection, removal of ties, etc.
    % Get all sequences
    seq_matrix = pt(whichPt).seq_matrix;
    
    % Remove sequences with too many ties
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

    % Get spikes in sequences
    all_spikes = [];
    all_times_all = [];
    for i = 1:size(seq_matrix,2)
        nonan = find(~isnan(seq_matrix(:,i)));
        times = seq_matrix(nonan,i);

        % Resort by time
        [~,I] = sort(times);
        nonan = nonan(I);

        all_spikes = [all_spikes;nonan];
        all_times_all = [all_times_all;seq_matrix(nonan,i)];
    end
    
    % Remove ictal spikes
    t = find(any(all_times_all >= (szTimes(:,1)-repmat(60,size(szTimes,1),1))' ...
        & all_times_all <= szTimes(:,2)',2));
    all_spikes(t) = [];
    fprintf('Removed %d ictal spikes \n',length(t));
    
    spike_in_seq_rate = length(all_spikes)/all_times/nchs*60;
    post_proc = [post_proc;spike_in_seq_rate];
    
    %% Post-clustering
    all_spikes = cluster(whichPt).all_spikes; % all spike channels
    idx = cluster(whichPt).idx;
    bad_cluster = cluster(whichPt).bad_cluster;
    bad_idx = find(ismember(idx,bad_cluster));
    all_spikes(bad_idx) = [];
    
    spikes_post_cluster = length(all_spikes)/all_times/nchs*60;
    post_clust = [post_clust;spikes_post_cluster];
    
    %fprintf('%s had %d spikes post cluster.\n',pt(whichPt).name,length(all_spikes));
    all_spike_num = [all_spike_num;length(all_spikes)];
    
end

fprintf('Before processing, the spike rate in spikes/ch/min was:\n%1.3f (range %1.3f-%1.3f)\n',...
    mean(all),min(all),max(all));

fprintf('After processing, the spike rate in spikes/ch/min was:\n%1.3f (range %1.3f-%1.3f)\n',...
    mean(post_proc),min(post_proc),max(post_proc));

fprintf('After clustering, the spike rate in spikes/ch/min was:\n%1.3f (range %1.3f-%1.3f)\n',...
    mean(post_clust),min(post_clust),max(post_clust));

fprintf('The average time analyzed was %1.1f hours (range %1.1f-%1.1f).\n',...
    mean(times_all_pts)/3600,min(times_all_pts)/3600,max(times_all_pts)/3600);

fprintf('After clustering, the mean number of spikes evaluated was:\n %d (range %d-%d).\n',...
    round(mean(all_spike_num)),min(all_spike_num),max(all_spike_num));

end