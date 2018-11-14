function chConnections(pt,whichPt)


%% Other ideas
%{
-Plot SOZ spike/seq freq over time
- plot seq frequence of leader chs over time
%}

%% Goal
%{
To see if the electrodes that have the biggest volume of influence (ie,
volume of electrodes that spike after the electrode of interest in a spike
sequence) are the seizure onset zone electrodes

I would expect that this would likely be the case for seizures that
generalize, but probably not for focal seizures
%}

%% Parameters
map_text = 'jet';
approach = 1;
doBootstrap = 0;
alpha = 99.9;
newSOZ =1;


%% Get sequences
%{

%}
all_times = pt(whichPt).cluster.all_times;
all_seq_cat = pt(whichPt).cluster.all_seq_cat;
bad_cluster = pt(whichPt).cluster.bad_cluster;
idx = pt(whichPt).cluster.idx;
bad_idx = find(ismember(idx,bad_cluster));
all_times(bad_idx) = [];
all_seq_cat(:,bad_idx) = [];

locs = pt(whichPt).electrodeData.locs(:,2:4);
nchs = size(locs,1);


%{
%% Get sequences and Remove sequences with too many ties
[all_seq_cat,all_times,~,~] = divideIntoSzChunksGen(pt,whichPt);
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
%}

%% Get colormap
fh_map = str2func(map_text);

%% Get SOZ channels
if newSOZ == 0
    soz = [];
    for j = 1:length(pt(whichPt).sz)
        soz = [soz;pt(whichPt).sz(j).chs];
    end
    soz = unique(soz);
else
   soz = pt(whichPt).newSOZChs; 
end

%% Get sz times
szTimes = zeros(length(pt(whichPt).sz),1);
for j = 1:length(pt(whichPt).sz)
   szTimes(j) = pt(whichPt).sz(j).onset;
end


%% Construct a matrix of channel connections
chCh = zeros(nchs,nchs);
for i = 1:size(all_seq_cat,2)
    seq = all_seq_cat(:,i);
    spike_chs = find(isnan(seq) == 0);
    spike_times = seq(find(isnan(seq) == 0));
    
    % sort channels by spike time
    [spike_times,I] = sort(spike_times);
    spike_chs = spike_chs(I);
    
    
    if approach == 1
        % Getting all unidirectional pairwise connections
        B = nchoosek(spike_chs,2);
        
        for j = 1:size(B,1)
        % I am saying that chCh(i,j) is the number of times that channel i
        % influences channel j
        chCh(B(j,1),B(j,2)) = chCh(B(j,1),B(j,2)) + 1;
        end
    elseif approach == 2
        % Get all unidirectional pairwise connections that start with first
        % channel
        B = nchoosek(spike_chs,2);
        B=B(B(:,1)==spike_chs(1),:);
        
        for j = 1:size(B,1)
        chCh(B(j,1),B(j,2)) = chCh(B(j,1),B(j,2)) + 1;
        end
    elseif approach == 3
        % Just get connection from first to 2nd channel
        chCh(spike_chs(1),spike_chs(2)) = ...
            chCh(spike_chs(1),spike_chs(2)) + 1;
        
    end
    
    
    
end


% Plot the pairwise connections
if 1 == 1
figure
imagesc(chCh)
colorbar
end

%% Get sequence number per channel
seq_freq = nansum(all_seq_cat,2);
leader_freq = sum(chCh,2);

%% Get significant connections
if doBootstrap == 1
    % Bootstrap approach
    ncons = sum(sum(chCh));
    fprintf('There are %d total connections.\n',ncons);
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

    if 1 == 0
        figure
        imagesc(squeeze(mean(chCh_all,1)))
    end

    s_con = sort(chCh_all(:));
    %scatter(1:length(s_con),s_con)
    perc = prctile(s_con,alpha);
    minCount = perc;
    fprintf(['By permutation testing, the minimum number of counts for a\n'...
        'connection to be significant is\n'...
        '%d for an alpha of %1.1fth percentile\n\n'],perc,alpha);
    
    ncons = sum(sum(chCh));
    lambda = ncons/nchs^2;
    X = poissinv(alpha/100,lambda);
    fprintf(['By poisson assumption, the number of counts is:\n'...
        '%d\n\n'],X);
    
    
else
    % Assume poisson distribution
    
    ncons = sum(sum(chCh));
    lambda = ncons/nchs^2;
    X = poissinv(alpha/100,lambda);
    minCount = X;
    
    fprintf(['Not doing permutation test, assuming poisson distribution.\n'...
        'Doing so yields min count number for significance of\n'...
        '%d for an alpha of %1.1fth percentile\n\n'],X,alpha);
end


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

gs = (fh_map(50));
[Y,E] = discretize(log(chull),size(gs,1));
figure
subplot(2,2,1)
set(gcf,'Position',[200 200 1000 700]);
scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
hold on

scatter3(locs(isnan(Y)==0,1),locs(isnan(Y)==0,2),locs(isnan(Y)==0,3),100,gs(Y(isnan(Y)==0),:));
[~,I] = max(chull);
scatter3(locs(I,1),locs(I,2),locs(I,3),100,gs(Y(I),:),'filled');
scatter3(locs(soz,1),locs(soz,2),locs(soz,3),'*','k');
title('Volume of convex hull of influenced channels');


all_con = chCh(:);
con_col = log(all_con);
gs = fh_map(round(max(con_col)));
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

gs = fh_map(50);
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

gs = fh_map(50);
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


%% Plot SOZ seq freq and max hull seq freq over time

% Get time range
all_times = floor(min(all_seq_cat(:,1))):ceil(max(all_seq_cat(:,end)));

% Get spikes containing SOZ
if length(soz) > 1
    seq_w_soz = any(isnan((all_seq_cat(soz,:)))==0,1);
else
    seq_w_soz = (isnan((all_seq_cat(soz,:)))==0);
end
times_w_soz = min(all_seq_cat(:,seq_w_soz));

% Get spikes starting in the SOZ
[~,firstCh] = min(all_seq_cat,[],1);
seq_s_soz = (ismember(firstCh,soz));
times_s_soz = min(all_seq_cat(:,seq_s_soz));


% Get spikes starting in the hull
[~,I] = max(chull);
seq_s_hull = (firstCh==I);
times_s_hull = min(all_seq_cat(:,seq_s_hull));

% Get spikes containing the hull
seq_w_hull = (isnan(all_seq_cat(I,:))==0);
times_w_hull = min(all_seq_cat(:,seq_w_hull));


% Plot scatter
figure
s1=scatter(times_w_soz,4*ones(length(times_w_soz),1),100);
hold on
s2=scatter(times_s_soz,3*ones(length(times_s_soz),1),100,'r');
s3=scatter(times_w_hull,2*ones(length(times_w_hull),1),100,'g');
s4=scatter(times_s_hull,1*ones(length(times_s_hull),1),100,'m');
yl = ylim;
for j = 1:length(szTimes)
   plot([szTimes(j) szTimes(j)],yl,'k','LineWidth',2); 
end
legend([s1 s2 s3 s4],{'Contains SOZ','Started with SOZ','Contains max hull','Started max hull'});

% Plot moving average of counts
window = 3600;

figure


subplot(4,1,1)
[counts,window_times] = movingSumCounts(times_w_soz,all_times,window);
plot(window_times,counts);
hold on
yl = ylim;
for j = 1:length(szTimes)
   plot([szTimes(j) szTimes(j)],yl,'k','LineWidth',2); 
end
title('Number of sequences containing SOZ');

subplot(4,1,2)
[counts,window_times] = movingSumCounts(times_s_soz,all_times,window);
plot(window_times,counts);
hold on
yl = ylim;
for j = 1:length(szTimes)
   plot([szTimes(j) szTimes(j)],yl,'k','LineWidth',2); 
end
title('Number of sequences starting in SOZ')

subplot(4,1,3)
[counts,window_times] = movingSumCounts(times_w_hull,all_times,window);
plot(window_times,counts);
hold on
yl = ylim;
for j = 1:length(szTimes)
   plot([szTimes(j) szTimes(j)],yl,'k','LineWidth',2); 
end

title('Number of sequences containing max connected electrode');
subplot(4,1,4)
[counts,window_times] = movingSumCounts(times_s_hull,all_times,window);
plot(window_times,counts);
hold on
yl = ylim;
for j = 1:length(szTimes)
   plot([szTimes(j) szTimes(j)],yl,'k','LineWidth',2); 
end
title('Number of sequences starting in max connected electrode')


end