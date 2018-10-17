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
   seq_tracker(isnan(allSeqs)== 0) = 1;
   num_seqs = nansum(seq_tracker,2);
   perc_seqs = num_seqs/nansum(num_seqs);
   [S,I] = sort(perc_seqs,'descend');
   CS = cumsum(S);
   
    
end


end

