function leadInteric(pt,whichPts)

for whichPt = whichPts
    
seq_matrix = pt(whichPt).seq_matrix;    
locs = pt(whichPt).electrodeData.locs(:,2:4);
soz = pt(whichPt).newSOZChs;

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
lambda = sum(sum(~isnan(seq_matrix)))/size(seq_matrix,1);
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
gs = (gray(50));
figure
scatter3(locs(:,1),locs(:,2),locs(:,3),60,'k');
[Ylat,~] = discretize(lat_sig,size(gs,1));
hold on
scatter3(locs(sig_ch,1),locs(sig_ch,2),locs(sig_ch,3),...
    60,gs(Ylat,:),'filled');
scatter3(locs(soz,1),locs(soz,2),locs(soz,3),200,'r+');
scatter3(locs(early_sig,1),locs(early_sig,2),locs(early_sig,3),...
    100,'g+');
    
end



end