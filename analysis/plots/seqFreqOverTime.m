function [chunk_seqs,times_plot,MI,rl] = seqFreqOverTime(pt,whichPt,window)

%% Parameters
dmin = 31;
nchs = length(pt(whichPt).channels);

% output file name
[~,~,~,resultsFolder,~] = fileLocations;

%% Get wij
xyChan = pt(whichPt).electrodeData.locs;
wij = getwij(xyChan,dmin);

%% First need to divide it up into multiple seizure chunks for plotting

% just look at first seizure
sequences = pt(whichPt).sz(1).seq_matrix;
sequences(sequences==0) = nan; % WHY ARE THERE ANY ZEROS?????


% Fill up the info from the first seizure
szTime = pt(whichPt).sz(1).onset;
seq_all{1} = sequences;

% Loop through the other seizures
for j = 2:length(pt(whichPt).sz)
    szTimeLast = pt(whichPt).sz(j-1).onset;
    szTime = pt(whichPt).sz(j).onset;
    totalTime = pt(whichPt).sz(j).runTimes(end,2) - pt(whichPt).sz(j).runTimes(1,1);
    
    sequences = pt(whichPt).sz(j).seq_matrix;
    sequences(sequences==0) = nan; % WHY ARE THERE ANY ZEROS?????
    
    if szTime - szTimeLast > totalTime
        % If the seizure time is more than 24 hours after the last seizure,
        % put this in a new chunk for plotting purposes
        seq_all{end+1} = sequences;
    else
        % if it is less than 12 hours after the last seizure, then need to
        % add any spikes that occured after the last seizure run time to
        % this new chunk
        keepAfter = pt(whichPt).sz(j-1).runTimes(end,2);
        firstSpikes = min(sequences,[],1);  
        seq_to_keep = sequences(:,firstSpikes > keepAfter);
        seq_all{end} = [seq_all{end},seq_to_keep];
    end

end

chunk_seqs = cell(size(seq_all));
chunk_seqs_chs = cell(size(seq_all));
times_plot = cell(size(seq_all));
rl = cell(size(seq_all));
MI = cell(size(seq_all));

%% Now divide the sequences into windows
for i = 1:length(seq_all)
   seq = seq_all{i};
   firstSpikes = min(seq,[],1);
   totalTime = firstSpikes(end) - firstSpikes(1);
   nchunks = ceil(totalTime/window);
   
   chunk_seqs{i} = zeros(nchunks,1);
   chunk_seqs_chs{i} = zeros(nchunks,nchs);
   times_plot{i} = zeros(nchunks,1);
   rl{i} = zeros(nchunks,nchs);
   MI{i} = zeros(nchunks,1);
   
   
   for tt = 1:nchunks
      times =  [(tt-1)*window + firstSpikes(1),tt*window + firstSpikes(1)];
      times_plot{i}(tt) = (times(1)+times(2))/2;
      
      % Get the appropriate sequences in this time
      correct_seqs = seq(:,firstSpikes >= times(1) & firstSpikes <= times(2));
      chunk_seqs{i}(tt) = size(correct_seqs,2);
      
      % get the number of sequences per channel (the starting channel???)
      for k = 1:size(correct_seqs,2)
         [~,ch] =  min(correct_seqs(:,k));
         chunk_seqs_chs{i}(tt,ch) = chunk_seqs_chs{i}(tt,ch) + 1;
          
      end
      
      %% Get recruitment latency for these sequences
      
      % for each sequence, the latency with which each channel is activated
      % in the sequence is the spike time in that channel minus the spike
      % time in the channel activated the earliest in that sequence
      latency_all_seq = correct_seqs - min(correct_seqs,[],1);
      
      % Take the average latency for the channel over all sequences
      mean_latency = nanmean(latency_all_seq,2);
      rl{i}(tt,:) = mean_latency';
      
      % Get the moran index
      MIstruct= moranStats(mean_latency',wij,nchs);
      if MIstruct.I > 1
          error('look\n');
      end
      MI{i}(tt) = MIstruct.I;
      
   end
   
   chunk_seqs{i} = chunk_seqs{i}/window;
   chunk_seqs_chs{i} = chunk_seqs_chs{i}/window;
   
    
end





%{
figure
for i = 1:length(chunk_seqs)
    times = times_plot{i};
    plot(times/3600,chunk_seqs{i},'k');
    hold on
end

yl = ylim;

for j = 1:length(pt(whichPt).sz)
    szTimes = [pt(whichPt).sz(j).onset,pt(whichPt).sz(j).offset];
    meanSzTimes = (szTimes(1) + szTimes(2))/2;
    plot([meanSzTimes meanSzTimes]/3600,[yl(1) yl(2)],'k--');
end

xlabel('Hour');
ylabel('Sequences per hour');
title('Sequence frequency over time');
set(gca,'FontSize',15)
%}

end