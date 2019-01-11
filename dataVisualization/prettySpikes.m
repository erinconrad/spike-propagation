function prettySpikes(pt,cluster)




%% Parameters
dur = [-1 1];

% HUP078
offset = [-3 27 7.7653];
%offset = [-3 40.5780 7.7653];
%offset = [-8.5508 40.5780 7.7653]; CORRECT

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

outputFolder = [resultsFolder,'pretty_plots/Fig1/'];

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

%% Initialize figure
figure
[ha,pos] = tight_subplot(1,3,[0.01 0.01],[0.01 0.01],[0.01 0.01]);

%% Plot EEG tracings
axes(ha(1));
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
axis(ha(1),'off');

%% Plot brain
axes(ha(2))
locs = pt(whichPt).electrodeData.locs(:,2:4);

% Get transformation matrix to get new coordinate locations
A = makeNewElecData(pt,whichPt);
%offset = [-10 30 0]; % bs
newlocs = A*locs-offset;

% Load gifti
brainFolder = '/Users/erinconrad/Desktop/residency stuff/R25/actual work/data/brains/';
giftiFolder = [brainFolder,pt(whichPt).name,'/'];
names = dir([giftiFolder,'*pial.gii']);
fname2 = names(1).name;
g = gifti([giftiFolder,fname2]);

circSize = 130;
p = plotGIFTI(g);
hold on
view(-120,-11);
%alpha(p,0.4)
scatter3(newlocs(:,1),newlocs(:,2),newlocs(:,3),circSize,'k','filled');
hold on

for i = 1:length(times)
    scatter3(newlocs(chs(i),1),newlocs(chs(i),2),newlocs(chs(i),3),circSize,...
        c_idx(i,:),'filled');   
    %{
   text(newlocs(chs(i),1),...
        newlocs(chs(i),2),...
        newlocs(chs(i),3),...
        sprintf('Spike %d',i),'FontSize',25,'color',c_idx(i,:));
        %}
   
   
end
xticklabels([]);
yticklabels([]);
zticklabels([]);

%view(-112,-11);


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
axes(ha(3))
yPos = [0 100 200];
groupRange = zeros(3,2);
textPos = {'X','Y','Z'};

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
    text(-0.7,(yLower+ywidth*(j-1) + yLower+ywidth*(j))/2,...
        textPos{j},'color','k','FontSize',30);
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
axis(ha(3),'off');

%% Adjust figure
set(gcf,'Position',[71 241 1237 538]);

%% Plot arrows and text
annotation('textarrow',[0.1 0.555],[0.8 0.355],'String','',...
    'color','b','FontSize',20,'LineWidth',2);
annotation('textarrow',[0.13 0.57],[0.6 0.385],'String','',...
    'color','b','FontSize',20,'LineWidth',2);
annotation('textarrow',[0.155 0.52],[0.4 0.28],'String','',...
    'color','r','FontSize',20,'LineWidth',2);
annotation('textarrow',[0.19 0.475],[0.2 0.28],'String','',...
    'color','r','FontSize',20,'LineWidth',2);

annotation('textbox',[0.33 0.9 0.5 0.1],...
    'String','Spike location clustering and plotting',...
    'linestyle','none','FontSize',23);

annotation('textbox',[0.01 0.88 0.3 0.1],'String','A',...
    'linestyle','none','FontSize',30);

annotation('textbox',[0.3 0.88 0.3 0.1],'String','B',...
    'linestyle','none','FontSize',30);

annotation('textbox',[0.65 0.88 0.3 0.1],'String','C',...
    'linestyle','none','FontSize',30);

print(gcf,[outputFolder,'Fig1'],'-depsc');
eps2pdf([outputFolder,'Fig1','.eps']);

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