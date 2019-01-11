function prettySpikes(pt,cluster)


%% Parameters
dur = [-1 1];


%% Spike info
whichPt = 8;

spikes = [7622,7623,7640,7645];
times = cluster(whichPt).all_times_all(spikes);
chs = cluster(whichPt).all_spikes(spikes);
idx = cluster(whichPt).idx(spikes);
c_idx = zeros(length(spikes),3);
for i = 1:length(idx)
    if idx(i) == 1, c_idx(i,:) = [0 0 1]; elseif idx(i) == 2, c_idx(i,:) = [1 0 0]; end
end

%% Get stuff
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;

if exist('pretty_spikes.mat','file') ~= 0
    load('pretty_spikes.mat')
else

    fs = pt(whichPt).fs;

    for i = 1:length(times)
        desiredTimes = times(i) + dur;
        desiredIndices = round(desiredTimes(1)*fs):round(desiredTimes(2)*fs);
        data = getiEEGData(pt(whichPt).ieeg_name,pt(8).channels,desiredIndices,pwfile);
        plot_points{i} = data.values(:,chs(i));
    end

    save('pretty_spikes.mat','plot_points','spikes');
end

%% Plot EEG tracings
figure
toAdd = 0;
timeAdd = 0;
for i = 1:length(times)
   
     timePlot = linspace(0,2,length(plot_points{i})) + timeAdd;
     plot(timePlot,plot_points{i}+toAdd,'color',c_idx(i,:),'linewidth',2);
     hold on
     textSp = sprintf('Spike %d',i);
     text(timePlot(end)+0.02,plot_points{i}(end)+toAdd,textSp,...
         'fontsize',20,'color',c_idx(i,:));
     toAdd = toAdd + min(plot_points{i}) - max(plot_points{i});
     timeAdd = timeAdd + 0.3;
     
     
end
xlim([0 timePlot(end)+0.6]);
xticks([])
yticks([])

%% Plot brain
locs = pt(whichPt).electrodeData.locs(:,2:4);

% Get transformation matrix to get new coordinate locations
A = makeNewElecData(pt,whichPt);
offset = [-10 30 0]; % bs
newlocs = A*locs-offset;

% Load gifti
brainFolder = '/Users/erinconrad/Desktop/residency stuff/R25/actual work/data/brains/';
giftiFolder = [brainFolder,pt(whichPt).name,'/'];
names = dir([giftiFolder,'*pial.gii']);
fname2 = names(1).name;
g = gifti([giftiFolder,fname2]);

figure
p = plotGIFTI(g);
hold on
alpha(p,0.4)
scatter3(newlocs(:,1),newlocs(:,2),newlocs(:,3),100,'k','filled');
hold on
for i = 1:length(times)
    scatter3(newlocs(chs(i),1),newlocs(chs(i),2),newlocs(chs(i),3),100,...
        c_idx(i,:),'filled');   
end
xticklabels([]);
yticklabels([]);
zticklabels([]);


%{
% Get convex hull of channels in each cluster
color_clust = [0 0 1;1 0 0];
for i = 1:cluster(whichPt).k
    c_locs = cluster(whichPt).all_locs(cluster(whichPt).idx == i,:);
    DT = delaunayTriangulation(c_locs);
    [C,v] = convexHull(DT);
    
    
    a= trisurf(C,DT.Points(:,1),DT.Points(:,2),DT.Points(:,3), ...
       'FaceColor',color_clust(i,:));
   
    alpha(a,0.2);
end
%}

%% Plot x,y,z coordinates
% NEED TO FIX

yPos = [0 100 200];
groupRange = zeros(3,2);
textPos = {'X','Y','Z'};

figure
for j = 1:3
    timeAdd = 0;
    maxThing = 0;
    minThing = 2000;
    for i = 1:length(spikes)
        plot_thing = locs(chs(i),j) + yPos(j);
        if plot_thing > maxThing
            maxThing = plot_thing;
        end
        if plot_thing < minThing
            minThing = plot_thing;
        end
        scatter(timeAdd,plot_thing,100,c_idx(i,:),'filled');
        hold on
        timeAdd = timeAdd + 1;
    end
    groupRange(j,:) = [minThing maxThing];
    
end

ywidth = (groupRange(2,2)+groupRange(3,1))/2 - (groupRange(1,2)+groupRange(2,1))/2;
ylim([(groupRange(1,2)+groupRange(2,1))/2-ywidth,...
    (groupRange(2,2)+groupRange(3,1))/2+ywidth+10])

for j = 1:3
    yLower = get(gca,'Ylim'); yLower = yLower(1);
    text(-0.4,(yLower+ywidth*(j-1) + yLower+ywidth*(j))/2,...
        textPos{j},'color','k','FontSize',20);
end

for i = 1:length(spikes)
    h = text(i-1,groupRange(3,2) + 10,sprintf('Spike %d',i),...
        'color',c_idx(i,:),'FontSize',20);
    set(h,'Rotation',90);
end

xlim([-0.5 3.2]);

plot(get(gca,'xlim'),...
    [(groupRange(1,2)+groupRange(2,1))/2 (groupRange(1,2)+groupRange(2,1))/2],...
    'k--','Linewidth',2);

plot(get(gca,'xlim'),...
    [(groupRange(2,2)+groupRange(3,1))/2 (groupRange(2,2)+groupRange(3,1))/2],...
    'k--','Linewidth',2);

xticks([])
yticks([])


%{
figure

timeAdd = 0;
for i = 1:length(spikes)
    toAdd = 0;
    for j = 1:3
        plot_thing = locs(chs(i),j)+toAdd;
        scatter(timeAdd,plot_thing,100,c_idx(i,:),'filled');
        hold on
        toAdd = toAdd + 100;
    end
    timeAdd = timeAdd + 1;
end
%}


end