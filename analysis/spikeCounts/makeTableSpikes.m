function makeTableSpikes(pt,cluster)

%% Get patients to do it for
whichPts = [];
for i = 1:length(pt)
    if isempty(pt(i).seq_matrix) == 0
        if size(cluster(i).bad_cluster) < cluster(i).k
            whichPts = [whichPts,i];
        end
    end
end

% I should be doing all 20 of these patients
if isequal(whichPts,[1,4,6,7,8,9,12,14,15,16,17,18,19,20,22,24,25,27,30,31]) == 0
    error('Warning, not doing correct patients!\n');
end

%% Initialize variables
name = cell(length(whichPts),1);
duration_t = cell(length(whichPts),1);
spike_rate_t = cell(length(whichPts),1);
post_processing_t = cell(length(whichPts),1);
post_cluster_t = cell(length(whichPts),1);
all_spike_num_t = cell(length(whichPts),1);
sz_num_t = zeros(length(whichPts),1);

count = 0;
for whichPt = whichPts
    count = count + 1;
    
    
    %% Name
    name{count} = pt(whichPt).name;
    
    %% Basic info
    szTimes = pt(whichPt).newSzTimes;
    nchs = size(pt(whichPt).channels,1);
    
    %% Seizure number
    sz_num_t(count) = length(pt(whichPt).sz);
    
    %% Get spike rate before processing; also get duration
    % Get number of spikes
    n_spikes_all = pt(whichPt).stats.nspikes;
    
    % Get duration
    duration= sum(diff(pt(whichPt).allTimes,1,2));
    duration_t{count} = sprintf('%1.1f',duration/3600);
    
    % Spike rate
    spike_rate = n_spikes_all/duration/nchs*3600;
    spike_rate_t{count} = sprintf('%1.1f',spike_rate);
    
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
    
    spike_in_seq_rate = length(all_spikes)/duration/nchs*3600;
    post_processing_t{count} = sprintf('%1.1f',spike_in_seq_rate);
    
    %% Post-clustering
    all_spikes = cluster(whichPt).all_spikes; % all spike channels
    idx = cluster(whichPt).idx;
    bad_cluster = cluster(whichPt).bad_cluster;
    bad_idx = find(ismember(idx,bad_cluster));
    all_spikes(bad_idx) = [];
    
    post_cluster = length(all_spikes)/duration/nchs*3600;
    post_cluster_t{count} = sprintf('%1.1f',post_cluster);
    
    %fprintf('%s had %d spikes post cluster.\n',pt(whichPt).name,length(all_spikes));
    all_spike_num = length(all_spikes);
    all_spike_num_t{count} = sprintf('%1.1f',all_spike_num);
    
    
end

T = table(name,duration_t,spike_rate_t,post_processing_t,post_cluster_t,all_spike_num_t,...
    'VariableNames',{'ID','Duration','Raw_spike_rate',...
    'Post_processing_spike_rate','Post_clustering_spike_rate','Final_spike_number'})

%{
fprintf(['Average duration:  %1.1f (range %1.1f-%1.1f)\n'...
    'Average raw rate: %1.1f (range %1.1f-%1.1f)\n'...
    'Average post-processing rate:  %1.1f (range %1.1f-%1.1f)\n'...
    'Average post-clustering rate:  %1.1f (range %1.1f-%1.1f)\n'...
    'Average final number:  %1.1f (range %d-%d)\n'],...
    mean(
    %}
    
    sz_num_t
    
end