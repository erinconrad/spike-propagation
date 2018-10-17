function getSpikeDist(pt,whichPts)

for whichPt = whichPts
    
   allSeqs = [];
   
   %% Loop through seizures and gather all sequences
   for j = 1:length(pt(whichPt).sz)
      allSeqs = [allSeqs,pt(whichPt).sz(j).seq_matrix];     
   end
    
   %% Get numbers of sequences per channel
   nSeq = size(allSeqs,2);
   seq_tracker = allSeqs;
   
   % put a 1 in the position of any channel that is involved in the
   % sequence
   seq_tracker(isnan(allSeqs)== 0) = 1;
   
   % summing this across all sequences gets the number of sequences in
   % which that channel is involved
   num_seqs = nansum(seq_tracker,2);
   
   % This is the total percentage of sequence-spikes that this channel is
   % involved in
   perc_seqs = num_seqs/nansum(num_seqs);
   
   % sort in descending order to get the most involved channels
   [S,~] = sort(perc_seqs,'descend');
   
   % Do the cumulative sum to see what number of channels is needed to
   % reach an arbitary percentage of sequence-spikes
   CS = cumsum(S);
   
   % the lowest index for which the cumulative sum is >0.9 is the number of
   % channels needed to account for 90% of all sequence spikes 
   index_to_90 = min(find(CS > 0.9));
   
   % the percentage of channels needed to account for 90% of all sequence
   % spikes
   perc_to_90 = index_to_90/length(CS);
   
   fprintf(['For %s, %d channels (%1.2f percent of %d channels) account ',...
       'for 90 percent of all sequence spikes\n'], pt(whichPt).name,...
       index_to_90, perc_to_90, length(CS));
    
end


end

