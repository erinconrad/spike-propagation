function permVecDist(sequences,electrodeData)

nperm = 1e2;
alpha = 0.05;
nseq = size(sequences,2);
dist_sum = zeros(nperm,1);

[~,true_dist] = getAngles(sequences,electrodeData);

true_dist_sum = sum(true_dist);

for i = 1:nperm
    dist = getDistPermute(sequences,electrodeData);
    dist_sum(i) = sum(dist);
    
end

[~,dist_sum] = sort(dist_sum);
CI = [dist_sum(max(1,round(alpha*nperm/2))),...
    dist_sum(min(length(dist_sum),round(length(dist_sum)-alpha*nperm/2)))];


end


function dist = getDistPermute(sequences,electrodeData)

nseq = size(sequences,2);
dist = zeros(nseq,1);

for i = 1:nseq
   curr_seq = sequences(:,i); 
   [B,I] = sort(curr_seq);
   C = I(isnan(B) == 0);
   nspikes = length(C);
   
   p = randperm(length(C));
   C = C(p);
   
   % get locations of the spikes
   locs = electrodeData.locs(C,2:4);
   
   % if odd number of spikes, ignore the middle spike
   if mod(nspikes,2) == 1
       early_idx = 1:(nspikes-1)/2;
       late_idx = (nspikes-1)/2+2:nspikes;
   else
       early_idx = 1:nspikes/2;
       late_idx = nspikes/2+1:nspikes;
   end
    
   % get mean coordinates of early channels and late channels
   early_mean = mean(locs(early_idx,:),1);
   late_mean = mean(locs(late_idx,:),1);
   
   % get distance between early and late
   dist(i) = norm(late_mean - early_mean);
   
end

end