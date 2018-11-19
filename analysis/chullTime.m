function chullTime(pt,whichPts)
map_text = 'jet';
fh_map = str2func(map_text);

for whichPt = whichPts
all_times = pt(whichPt).cluster.all_times;
all_seq_cat = pt(whichPt).cluster.all_seq_cat;
bad_cluster = pt(whichPt).cluster.bad_cluster;
idx = pt(whichPt).cluster.idx;
bad_idx = find(ismember(idx,bad_cluster));
all_times(bad_idx) = [];
all_seq_cat(:,bad_idx) = [];
nseq = length(all_times);

locs = pt(whichPt).electrodeData.locs(:,2:4);
nchs = size(locs,1);


%% Get SOZ channels
soz = pt(whichPt).newSOZChs; 

%% Get sz times
szTimes = zeros(length(pt(whichPt).sz),1);
for j = 1:length(pt(whichPt).sz)
   szTimes(j) = pt(whichPt).sz(j).onset;
end


chull = zeros(nseq,1);
plength = zeros(nseq,1);
chull_ch = cell(nchs,1);
for i = 1:nseq
    seq = all_seq_cat(:,i);
    spike_chs = find(isnan(seq) == 0);
    spike_times = seq(find(isnan(seq) == 0));
    
    
    for j = 1:length(spike_chs)-1
        plength(i) = plength(i) + ...
            vecnorm(locs(spike_chs(j),:) - locs(spike_chs(j+1),:));
        
    end
    
    if length(spike_chs) <4, continue; end
    
    % sort channels by spike time
    [spike_times,I] = sort(spike_times);
    spike_chs = spike_chs(I);
    spike_locs = locs(spike_chs,:);
    [K,V] = convhull(spike_locs(:,1),spike_locs(:,2),spike_locs(:,3));
    chull(i) = V;
    chull_ch{spike_chs(1)} = [chull_ch{spike_chs(1)} V];
    
    
end

median_chull = zeros(nchs,1);

for i = 1:nchs
   median_chull(i) = median(chull_ch{i}); 
end
median_chull(isnan(median_chull)) = 0;

[~,I] = max(median_chull);
figure
gs = fh_map(50);
[Y,E] = discretize(median_chull,size(gs,1));
scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
hold on
scatter3(locs(median_chull~=0,1),...
    locs(median_chull~=0,2),...
    locs(median_chull~=0,3),100,gs(Y(median_chull~=0),:));
scatter3(locs(I,1),locs(I,2),locs(I,3),100,gs(Y(I),:),'filled');
scatter3(locs(soz,1),locs(soz,2),locs(soz,3),'*');


if 1 == 0
window = 600;
nbins = round((all_times(end) - all_times(1))/window);
[Y,E] = discretize(all_times,nbins);
bin_times = zeros(nbins,1);
bin_vol = zeros(nbins,1);
for i = 1:nbins
   if isempty(all_times(Y==i)) == 1
        bin_times(i) = nan;
        bin_vol(i) = nan;
        continue
   end
   bin_times(i) = max(all_times(Y==i)); 
   bin_vol(i) = median(chull(Y==i));
end

%scatter(1:nseq,chull)
%scatter(all_times/3600,(chull))
scatter(bin_times/3600,bin_vol)
%scatter(all_times/3600,plength)
%hold on
%plot(all_times/3600,smooth(plength,100))
hold on
%plot(all_times/3600,log(chull))
for j = 1:size(szTimes,1)
   yl = ylim; 
   szOnset = szTimes(j,1);
   sz = plot([szOnset szOnset]/3600,yl,'k','LineWidth',3);
end

end
    
end

end