function chConnections(pt,whichPt)


%% Goal
%{
To see if the electrodes that have the biggest volume of influence (ie,
volume of electrodes that spike after the electrode of interest in a spike
sequence) are the seizure onset zone electrodes

I would expect that this would likely be the case for seizures that
generalize, but probably not for focal seizures
%}

%% Parameters

alpha = 99.9;


%% Get sequences
[all_seq_cat,all_times,~,~] = divideIntoSzChunksGen(pt,whichPt);
locs = pt(whichPt).electrodeData.locs(:,2:4);
nchs = size(locs,1);


%% Remove sequences with too many ties
keep = ones(size(all_seq_cat,2),1);
for s = 1:size(all_seq_cat,2)
   curr_seq = all_seq_cat(:,s);
   nonans = curr_seq(~isnan(curr_seq));
   norepeats = unique(nonans);
   if length(norepeats) < 0.5*length(nonans)
       keep(s) = 0;
   end
end

all_seq_cat(:,keep==0) = [];
all_times(:,keep == 0) = [];

fprintf(['%s had %d sequences (%1.2f of all sequences) deleted'...
    'for having >50 percent ties\n%d sequences remain\n'],...
    pt(whichPt).name,sum(keep == 0),sum(keep == 0)/length(keep),sum(keep==1));



%% Get SOZ channels
soz = [];
for j = 1:length(pt(whichPt).sz)
    soz = [soz;pt(whichPt).sz(j).chs];
end
soz = unique(soz);

%% Construct a matrix of channel connections
chCh = zeros(nchs,nchs);
for i = 1:size(all_seq_cat,2)
    seq = all_seq_cat(:,i);
    spike_chs = find(isnan(seq) == 0);
    spike_times = seq(find(isnan(seq) == 0));
    
    % sort channels by spike time
    [spike_times,I] = sort(spike_times);
    spike_chs = spike_chs(I);
    
    % Note that doing this generates a DIRECTIONAL bunch of channel pairs.
    % It always be (1st ch, 2nd ch) and not (2nd ch, 1st ch)
    B = nchoosek(spike_chs,2);
    
    for j = 1:size(B,2)
        % I am saying that chCh(i,j) is the number of times that channel i
        % influences channel j
        chCh(B(j,1),B(j,2)) = chCh(B(j,1),B(j,2)) + 1;
    end
    
end


% Plot the pairwise connections
if 1 == 0
imagesc(chCh)
colorbar
end

%% Get sequence number per channel
seq_freq = nansum(all_seq_cat,2);
leader_freq = sum(chCh,2);

%% Bootstrap to get significance
if 1 == 1
ncons = sum(sum(chCh));
nboot = 1e3;
max_size = nchs*nchs;
chCh_all = zeros(nboot,nchs,nchs);
for ib = 1:nboot
    if mod(ib,100) == 0
        fprintf('Doing %d of %d\n', ib,nboot);
    end
    chCh_f = zeros(nchs,nchs);
    for j = 1:ncons
       chCh_f(randi(max_size)) = chCh_f(randi(max_size)) + 1;
    end
    chCh_all(ib,:,:) = chCh_f;
end


figure
imagesc(squeeze(mean(chCh_all,1)))
s_con = sort(chCh_all(:));
scatter(1:length(s_con),s_con)
perc = prctile(s_con,alpha);
end

minCount = perc;
fprintf(['By permutation testing, the minimum number of counts for a\n'...
    'connection to be significant is\n'...
    '%d for an alpha of %dth percentile\n'],perc,alpha);

%% Get convex hull of the influence of each channel

% For each channel, get list of channels that are influenced by it enough
chInfluence = cell(nchs,1);

for i = 1:length(chInfluence)
    for j = 1:size(chCh,2)
        if chCh(i,j) > minCount % get better threshold
            chInfluence{i} = [chInfluence{i},j];
        end
    end   
end

% Calculate convex hull volume 
chull = zeros(nchs,1);
for i = 1:nchs
   if isempty(chInfluence{i}) == 1, continue; end
   hull_chs = chInfluence{i};
   if length(hull_chs) <= 2, continue; end
   hull_locs = locs(hull_chs,:);
   hull_locs = [locs(i,:);hull_locs];
   [K,V] = convhull(hull_locs(:,1),hull_locs(:,2),hull_locs(:,3));
   
   chull(i) = V;
   %{
   figure
   set(gcf,'Position',[200 200 800 550]);
   scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
   hold on
   scatter3(locs(i,1),locs(i,2),locs(i,3),100,'r');
   plot3(hull_locs(K,1),hull_locs(K,2),hull_locs(K,3));
   %}
end

gs = parula(50);
[Y,E] = discretize(log(chull),size(gs,1));
figure
subplot(2,2,1)
set(gcf,'Position',[200 200 1000 700]);
scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
hold on

scatter3(locs(isnan(Y)==0,1),locs(isnan(Y)==0,2),locs(isnan(Y)==0,3),100,gs(Y(isnan(Y)==0)));
[~,I] = max(chull);
scatter3(locs(I,1),locs(I,2),locs(I,3),100,gs(Y(I)),'filled');
scatter3(locs(soz,1),locs(soz,2),locs(soz,3),'*','k');
title('Volume of convex hull of influenced channels');


all_con = chCh(:);
con_col = log(all_con);
gs = parula(max(con_col));
[Y,E] = discretize(con_col,size(gs,1));
Y = reshape(Y,[size(chCh,1),size(chCh,1)]);

subplot(2,2,2)
scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
hold on
for i = 1:size(chCh,1)
   for j = 1:size(chCh,2) 
       con = chCh(i,j);
       dp = locs(j,:) - locs(i,:);
       if con > minCount
           
           quiver3(locs(i,1),locs(i,2),locs(i,3),dp(1), dp(2), dp(3),...
               'color',gs(Y(i,j),:));
           %{
           
           quiver3(locs(i,1),locs(i,2),locs(i,3),dp(1), dp(2), dp(3),...
               'color',[0 0 0],'linewidth',0.5);
               %}
       end
   end
end
scatter3(locs(soz,1),locs(soz,2),locs(soz,3),'*','k');
title('Connections between channels');

%% Seq freq

gs = parula(50);
[Y,E] = discretize(log(seq_freq),size(gs,1));
[~,most_freq] = max(seq_freq);


subplot(2,2,3)
scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
hold on
scatter3(locs(seq_freq~=0,1),locs(seq_freq~=0,2),locs(seq_freq~=0,3),...
    100,gs(Y(seq_freq~=0),:));
scatter3(locs(most_freq,1),locs(most_freq,2),locs(most_freq,3),...
    100,gs(Y(most_freq),:),'filled');
scatter3(locs(soz,1),locs(soz,2),locs(soz,3),'*','k');
title('Sequence frequency');

%% Leader freq

gs = parula(50);
[Y,E] = discretize(log(leader_freq),size(gs,1));
[~,most_freq] = max(leader_freq);


subplot(2,2,4)
scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
hold on
scatter3(locs(leader_freq~=0,1),locs(leader_freq~=0,2),locs(leader_freq~=0,3),...
    100,gs(Y(leader_freq~=0),:));
scatter3(locs(most_freq,1),locs(most_freq,2),locs(most_freq,3),...
    100,gs(Y(most_freq),:),'filled');
scatter3(locs(soz,1),locs(soz,2),locs(soz,3),'*','k');
title('Leader frequency');


%% Stats

% Distance from max convex hull to  sz onset mean
[~,I] = max(chull);
hullLoc = (locs(I,:));
szMeanLoc = mean(locs(soz,:),1);
hullDist = sqrt(sum((hullLoc-szMeanLoc).^2));

% Distance from all electrodes to sz onset mean
allCh = locs;
allChDist = sqrt(sum((allCh - repmat(szMeanLoc,nchs,1)).^2,2));

better = sum((allChDist<hullDist));
betterPerc = better/nchs*100;
fprintf('%1.1f percent of electrodes (%d of %d) outperformed our electrode\n',...
    betterPerc,better,nchs);



end