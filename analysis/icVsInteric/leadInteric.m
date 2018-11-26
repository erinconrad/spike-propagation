function leadInteric(pt,whichPts)

for whichPt = whichPts
    
seq_matrix = pt(whichPt).seq_matrix;    
locs = pt(whichPt).electrodeData.locs(:,2:4);

%% Remove sequences with too many ties
if 1 == 1
keep = ones(size(seq_matrix,2),1);
old_seq_matrix = seq_matrix;
for s = 1:size(seq_matrix,2)
   curr_seq = seq_matrix(:,s); 
   curr_seq = curr_seq(~isnan(curr_seq));
   norepeats = unique(curr_seq);
   if length(norepeats) < 0.5*length(curr_seq)
       keep(s) = 0;
   end
end

seq_matrix(:,keep==0) = [];
fprintf('Deleted %d (%1.1f of all sequences) for containing >50 percent ties\n',...
    sum(keep==0),sum(keep==0)/length(keep));

end


%% Get recruitment latency
lat = nanmean(seq_matrix - min(seq_matrix,[],1),2);

%% Get significant number of spikes and thus significant chs
lambda = sum(sum(~isnan(seq_matrix)))/size(seqs,1);
alpha1 = 99;
X = poissinv(alpha1/100,lambda);

% Get "significant" channels
nseq = sum(~isnan(seq_matrix),2);
sig_ch = find(nseq>X);

%% Get recruitment latency of significant channels
lat_sig = lat(sig_ch);
[~,earliest] = min(lat_sig);
early_sig = sig_ch(earliest);

%% Plot
figure
scatter3(
[Ylat,~] = discretize(lat_sig,
    
end



end