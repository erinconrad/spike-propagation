function [unclean,idx] = findUnclean(P,pt,sz,block)

nseq = size(P(pt).sz(sz).block(block).data.sequences,2)/2;
nclean = size(P(pt).sz(sz).block(block).data.cleanseq,2)/2;

full = P(pt).sz(sz).block(block).data.sequences;
full_clean = P(pt).sz(sz).block(block).data.cleanseq;

% Pad clean seqs with zeros
vertdiff = size(full,1)-size(full_clean,1);
full_clean = [full_clean;zeros(vertdiff,size(full_clean,2))];


unclean = [];
unclean_idx = 0;
idx = [];
for i = 1:nclean
   
   
   
   while 1
       
   temp_clean = full_clean(:,(i-1)*2+1:(i-1)*2+2);
   
   unclean_idx = unclean_idx + 1;
   temp_full = full(:,(unclean_idx-1)*2+1:(unclean_idx-1)*2+2);
       
   if isequal(temp_clean,temp_full) == 1
       break
   else
       idx = [idx,unclean_idx];
       unclean = [unclean,temp_full];
       
   end
   end
    
end


end