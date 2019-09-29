function aes_ex_seq_nogif(pt)


% This is for HUP078
whichPt = 8;
offset = [-3 27 7.7653];
whichSeq = 3200;
mark_size = 600;

%% Get stuff
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile,other] = fileLocations;
addpath(other.gifti)
addpath(genpath(scriptFolder))
outputFolder = [resultsFolder,'pretty_plots/aes_gif/'];
mkdir(outputFolder)
other_file_out = [outputFolder,'seq_aes_still'];


locs = pt(whichPt).electrodeData.locs(:,2:4);
seq = pt(whichPt).seq_matrix(:,whichSeq);

% Get the non nan times for the sequence
spike_chs = find(isnan(seq) == 0);
spike_times = seq(find(isnan(seq) == 0));
[spike_times,I] = sort(spike_times);
spike_chs = spike_chs(I);


% Get transformation matrix to get new coordinate locations
A = makeNewElecData(pt,whichPt);
%offset = [-10 30 0]; % bs
locs = A*locs-offset;
soz = pt(whichPt).newSOZChs;
szTimes = pt(whichPt).newSzTimes;
chs = 1:size(locs,1);

%% Load gifti
brainFolder = '/Users/erinconrad/Desktop/residency stuff/R25/actual work/data/brains/';
giftiFolder = [brainFolder,pt(whichPt).name,'/'];
names = dir([giftiFolder,'*pial.gii']);
fname2 = names(1).name;
g = gifti([giftiFolder,fname2]);


fig = figure;
set(fig,'position',[10 10 1000 800])
set(gcf,'color','white');

p = plotGIFTI(g);
hold on
%[0.6 0.6 0.6]
scatter3(locs(:,1),locs(:,2),locs(:,3),mark_size,'markerfacecolor',[0.3 0.3 0.3]);
scatter3(locs(spike_chs,1),locs(spike_chs,2),locs(spike_chs,3),mark_size,'markerfacecolor',...
    'w');
scatter3(locs(spike_chs(1),1),locs(spike_chs(1),2),locs(spike_chs(1),3),mark_size,'markerfacecolor',...
    [0 1 0]);
scatter3(locs(spike_chs(end),1),locs(spike_chs(end),2),locs(spike_chs(end),3),mark_size,'markerfacecolor',...
    [1 0 0]);
scatter3(locs(:,1),locs(:,2),locs(:,3),mark_size,'k','LineWidth',3); 



% Plot paths with arrows
for i = 1:length(spike_times)-1
    start = locs(spike_chs(i),:);
    finish = locs(spike_chs(i+1),:);
    dp = finish-start;
    quiver3([start(1)],[start(2)], [start(3)],...
        dp(1), dp(2), dp(3),'Color',[0 0 0],'LineWidth',4,'MaxHeadSize',4);    
end

if whichPt == 8
    view(-120,-11);
end

%print(fig,other_file_out,'-depsc');
%print(fig,other_file_out,'-dpng');



end