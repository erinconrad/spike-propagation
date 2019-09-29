function aes_aoi(pt,cluster)

whichPt = 8;

if whichPt == 8
    offset = [-3 27 7.7653];
elseif whichPt == 3
    offset = [-1.0029,2.3087,28.9465];
end
plotConn = 0;
removeTies = 1;
doBootstrap = 0;
alpha1 = 95;
map_text = 'jet';
fh_map = str2func(map_text);
new_bootstrap = 0;
fontsizes = 20;

[~,~,scriptFolder,resultsFolder,~,other] = fileLocations;
p1 = genpath(scriptFolder);
addpath(other.gifti)
addpath(p1);

destFolder = [resultsFolder,'influence/'];
mkdir(destFolder)

%% Patient parameters
fprintf('Doing %s\n',pt(whichPt).name);
locs = pt(whichPt).electrodeData.locs(:,2:4);
nchs = size(locs,1);
szTimes = pt(whichPt).newSzTimes;
soz = pt(whichPt).newSOZChs; 
saveFolder = [destFolder,pt(whichPt).name,'/'];
mkdir(saveFolder);

seq_matrix = pt(whichPt).seq_matrix;


%% Remove ties
if removeTies == 1
    keep = ones(size(seq_matrix,2),1);
    for s = 1:size(seq_matrix,2)
       curr_seq = seq_matrix(:,s);
       nonans = curr_seq(~isnan(curr_seq));
       norepeats = unique(nonans);
       if length(norepeats) < 0.5*length(nonans)
           keep(s) = 0;
       end
    end
    seq_matrix(:,keep==0) = [];
    fprintf(['%s had %d sequences (%1.2f of all sequences) deleted'...
    'for having >50 percent ties\n%d sequences remain\n'],...
    pt(whichPt).name,sum(keep == 0),sum(keep == 0)/length(keep),sum(keep==1));

end

%% Remove ictal sequences
all_times = seq_matrix(:);
icTimes = find(any(all_times >= (szTimes(:,1)-repmat(60,size(szTimes,1),1))' ...
    & all_times <= szTimes(:,2)',2));
seq_matrix(icTimes) = nan;
fprintf('Removed %d ictal spikes\n',length(icTimes));
%{
first_time = min(seq_matrix,[],1);
t = (any(first_time >= (szTimes(:,1)-repmat(60,size(szTimes,1),1)) ...
    & first_time <= szTimes(:,2),2));
seq_matrix(:,t) = [];
fprintf('Removed %d ictal spikes \n',sum(t));
%}

%% Get cluster info
all_times_all = cluster(whichPt).all_times_all; % all spike times
all_spikes = cluster(whichPt).all_spikes; % all spike channels
all_locs = cluster(whichPt).all_locs;
k = cluster(whichPt).k; % the number of clusters
idx = cluster(whichPt).idx; % the cluster index for every spike
C = cluster(whichPt).C; % the centroids of the clusters
bad_cluster = cluster(whichPt).bad_cluster; % which clusters are bad


%% Compare number of spikes in cluster array and my data
if sum(sum(~isnan(seq_matrix))) ~= length(all_times_all)
    error('Warning, number of spikes do not align\n');
end



%% Find bad spikes
bad_idx = find(ismember(idx,bad_cluster));

% Nx2 array of bad spikes, showing the channel and time
bad_spikes = [all_spikes(bad_idx),all_times_all(bad_idx)];

%% Get all sequences

new_seq_matrix = seq_matrix;
n_removed = 0;

%% Go through sequence matrix and remove bad spikes
for ich = 1:size(seq_matrix,1)
    % loop across electrodes

    % All spike times for this channel
    spikeTimesCh = seq_matrix(ich,:);

    % Get the bad spikes in that channel
    bad_times_for_ch = bad_spikes(bad_spikes(:,1) == ich,2);

    % Make sure I am finding all of them
    Lia = ismember(spikeTimesCh,bad_times_for_ch);
    if sum(Lia) ~= length(bad_times_for_ch)
        error(sprintf('Did not find all bad spikes for channel %d\n',ich));
    end

    %{
    if sum(Lia) > 0
        fprintf('Removed %d spikes for channel %d\n',sum(Lia),ich)
    end
    %}
    n_removed = n_removed + sum(Lia);

    % Make bad spikes nans
    spikeTimesCh(Lia==1) = nan;
    new_seq_matrix(ich,:) = spikeTimesCh;


end

if n_removed~=length(bad_idx)
    error('Incorrect number of bad spikes removed\n');
end
fprintf('Removed %d spikes for being in bad clusters\n',n_removed);


%% Remove sequences that have fewer than 5 spikes
removeSeq = zeros(size(new_seq_matrix,2),1);
for s = 1:size(new_seq_matrix,2)
    currSeq = new_seq_matrix(:,s);
    currSeq(isnan(currSeq)) = [];
    if length(currSeq) < 5
        removeSeq(s) = 1;
    end
end

fprintf('Removed %d sequences for now being too short\n',sum(removeSeq));
new_seq_matrix(:,removeSeq==1) = [];


seq_matrix = new_seq_matrix;

fprintf('%d sequences remain\n',size(seq_matrix,2));


%% Construct a matrix of channel connections
chCh = zeros(nchs,nchs);
all_spike_times = [];
upDown = zeros(nchs,2);
for i = 1:size(seq_matrix,2)
    seq = seq_matrix(:,i);
    spike_chs = find(isnan(seq) == 0);
    spike_times = seq((isnan(seq) == 0));

    % sort channels by spike time
    [spike_times,I] = sort(spike_times);
    spike_chs = spike_chs(I);

    all_spike_times = [all_spike_times;spike_times];

    % Get all unidirectional pairwise connections that start with first
    % channel
    B = nchoosek(spike_chs,2);
    B = B(B(:,1) == spike_chs(1),:);

    for j = 1:size(B,1)
        chCh(B(j,1),B(j,2)) = chCh(B(j,1),B(j,2)) + 1;
    end

    % Add info for degree preference
    for j = 1:length(spike_chs)
        % Add the number of channels before it in the sequence
        upDown(spike_chs(j),1) = upDown(spike_chs(j),1) + j-1;

       upDown(spike_chs(j),2) = upDown(spike_chs(j),2) + length(spike_chs)-j;
    end

end


if whichPt == 8, exampleChCh = chCh; end



% Assume poisson distribution (produces same result as permutation
% test)
ncons = sum(sum(chCh));
lambda = ncons/nchs^2;
X = poissinv(alpha1/100,lambda);
minCount = X;
fprintf(['By poisson assumption, the number of counts is:\n'...
    '%d\n\n'],X);


if whichPt == 8, exampleX = X; end

%% Now find connections that are more frequent by chance
n_spikes = length(all_spike_times);
lambda_spikes = n_spikes/nchs;
X_spikes = poissinv(alpha1/100,lambda_spikes);
minCountSpikes = X_spikes;

n_spikes_ch = sum(~isnan(seq_matrix),2);
ch_w_spikes = find(n_spikes_ch>minCountSpikes);
spiker = n_spikes_ch>minCountSpikes;

%% Get the indices of the channels that are influenced by each other
chInfluence = cell(nchs,1);
for i = 1:length(chInfluence)
    for j = 1:size(chCh,2)
        if chCh(i,j) > minCount
            chInfluence{i} = [chInfluence{i},j];
        end
    end
end

%% Get the surface are of influence of each channel
sa = zeros(nchs,1);
for i = 1:nchs
    if isempty(chInfluence{i}) == 1, continue; end
    downstream = chInfluence{i};
    if length(downstream) <= 2, continue; end
    down_locs = locs(downstream,:);
    down_locs = [locs(i,:);down_locs];
    tri = delaunay(down_locs(:,1),down_locs(:,2));
    P = [down_locs(:,1),down_locs(:,2),down_locs(:,3)];
    v1 = P(tri(:,2),:)-P(tri(:,1),:);
    v2 = P(tri(:,3),:)-P(tri(:,2),:);
    cp = 0.5*cross(v1,v2);
    sa(i) = sum(sqrt(dot(cp,cp,2)));
    %sa(i) = areaConnecting(down_locs,0);
end

%% seq freq
seq_freq = sum(~isnan(seq_matrix),2);

%% How far is the channel with max area of influence from SOZ?

if whichPt == 8, exampleSA = sa; end
if whichPt == 8, exampleChInfluence = chInfluence; end


%% Load gifti
brainFolder = '/Users/erinconrad/Desktop/residency stuff/R25/actual work/data/brains/';
giftiFolder = [brainFolder,pt(whichPt).name,'/'];
names = dir([giftiFolder,'*pial.gii']);
fname2 = names(1).name;
g = gifti([giftiFolder,fname2]);

locs = pt(whichPt).electrodeData.locs(:,2:4);
locs = pt(whichPt).electrodeData.locs(:,2:4);

% Get transformation matrix to get new coordinate locations
A = makeNewElecData(pt,whichPt);
%offset = [-10 30 0]; % bs
locs = A*locs-offset;
soz = pt(whichPt).newSOZChs;
szTimes = pt(whichPt).newSzTimes;
chs = 1:size(locs,1);

%[~,I] = max(exampleSA);
I = 83;%69;


%% Figures
if  0
    downstream = exampleChInfluence{I};
figure
set(gcf,'position',[10 10 900 900])
circSize = 600;

set(gcf,'color','white');
    
p = plotGIFTI(g);
hold on
scatter3(locs(:,1),locs(:,2),locs(:,3),circSize,'markerfacecolor',[0.3 0.3 0.3]);
scatter3(locs(downstream,1),locs(downstream,2),locs(downstream,3),circSize,...
    'markerfacecolor','w');
scatter3(locs(:,1),locs(:,2),locs(:,3),circSize,'k','linewidth',2);
hold on
scatter3(locs(I,1),locs(I,2),locs(I,3),circSize,'g','filled');


for j = 1:length(downstream)
    dp = locs(downstream(j),:) - locs(I,:);
    quiver3(locs(I,1),locs(I,2),locs(I,3),dp(1), dp(2), dp(3),...
        'color','k','linewidth',3,'maxheadsize',0.4);
end
scatter3(locs(I,1),locs(I,2),locs(I,3),circSize,'g','filled');
xticklabels([])
yticklabels([])
zticklabels([])
title('Downstream electrodes')
set(gca,'FontSize',30)
if whichPt == 8
    view(-120,-11);
end
end

%{
figure
set(gcf,'position',[10 10 900 900])
set(gcf,'color','white');
circSize = 600;    
p = plotGIFTI(g);
hold on
downstream = [I,downstream];
down_locs = locs(downstream,:);
tri = delaunay(down_locs(:,1),down_locs(:,2));
cv = trisurf(tri,down_locs(:,1),down_locs(:,2),...
    down_locs(:,3),'facecolor','r');
alpha(cv,0.2) 
hold on
scatter3(locs(:,1),locs(:,2),locs(:,3),circSize,'markerfacecolor',[0.3 0.3 0.3]);
scatter3(locs(downstream,1),locs(downstream,2),locs(downstream,3),circSize,...
    'markerfacecolor','w');
scatter3(locs(:,1),locs(:,2),locs(:,3),circSize,'k','linewidth',2);
hold on
scatter3(locs(I,1),locs(I,2),locs(I,3),circSize,'g','filled');
xticklabels([])
yticklabels([])
zticklabels([])
title(sprintf('Area of influence for electrode %d',I));
if whichPt == 8
    view(-120,-11);
end
set(gca,'FontSize',fontsizes)
%}

if 0
figure
set(gcf,'position',[10 10 900 900])
circSize = 600;

set(gcf,'color','white');
    
p = plotGIFTI(g);
hold on
scatter3(locs(:,1),locs(:,2),locs(:,3),circSize,exampleSA,'filled')
scatter3(locs(:,1),locs(:,2),locs(:,3),circSize,'k','linewidth',2);
xticklabels([])
yticklabels([])
zticklabels([])
if whichPt == 8
    view(-120,-11);
end
colorbar('ticks',[min(exampleSA),max(exampleSA)],...
    'ticklabels',{'Low','High'})

title(sprintf('Area of influence'));
set(gca,'FontSize',30)
end


% Seq freq
if 1
figure
set(gcf,'position',[10 10 900 900])
[ha, ~] = tight_subplot(1, 1, [.1 .01],[.05 .02],[.03 .03]);
circSize = 600;

set(gcf,'color','white');
axes(ha(1))
p = plotGIFTI(g);
hold on
scatter3(locs(:,1),locs(:,2),locs(:,3),circSize,seq_freq,'filled')
scatter3(locs(:,1),locs(:,2),locs(:,3),circSize,'k','linewidth',2);
scatter3(locs(soz,1),locs(soz,2),locs(soz,3),circSize,'p',...
    'markerfacecolor','w','markeredgecolor','k','linewidth',2);
xticklabels([])
yticklabels([])
zticklabels([])
if whichPt == 8
    view(-120,-11);
end
colorbar('ticks',[min(seq_freq),max(seq_freq)],...
    'ticklabels',{'Low','High'},'location','east')

t = text(-10,0,90,'Spike rate','fontsize',40);
set(gca,'FontSize',30)
alpha(p,0.2)
end

% SOZ
if 0
figure
set(gcf,'position',[10 10 900 900])
circSize = 600;

set(gcf,'color','white');
    
p = plotGIFTI(g);
hold on
scatter3(locs(:,1),locs(:,2),locs(:,3),circSize,'markeredgecolor',[0.3 0.3 0.3],...
    'markerfacecolor',[0.3 0.3 0.3])
hold on
scatter3(locs(soz,1),locs(soz,2),locs(soz,3),circSize,'w','filled');
scatter3(locs(:,1),locs(:,2),locs(:,3),circSize,'k','linewidth',2);
xticklabels([])
yticklabels([])
zticklabels([])
if whichPt == 8
    view(-120,-11);
end
alpha(p,0.2)


title(sprintf('Seizure onset zone'));
set(gca,'FontSize',30)
end

    
end


