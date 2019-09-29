function aes_resec_elecs(pt)

% This is for HUP078
whichPt = 8;
offset = [-3 27 7.7653];
mark_size = 600;

%% Get stuff
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile,other] = fileLocations;
addpath(other.gifti)
addpath(genpath(scriptFolder))
outputFolder = [resultsFolder,'pretty_plots/aes_gif/'];
mkdir(outputFolder)
other_file_out = [outputFolder,'resec_elecs'];

% Get list of resected electrodes
resec_elecs = pt(whichPt).resecElecs;

locs = pt(whichPt).electrodeData.locs(:,2:4);
% Get transformation matrix to get new coordinate locations
A = makeNewElecData(pt,whichPt);
%offset = [-10 30 0]; % bs
locs = A*locs-offset;


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
scatter3(locs(:,1),locs(:,2),locs(:,3),mark_size,'markerfacecolor',[0.3 0.3 0.3]);
scatter3(locs(resec_elecs,1),locs(resec_elecs,2),locs(resec_elecs,3),mark_size,'markerfacecolor',...
    'w');
scatter3(locs(resec_elecs,1),locs(resec_elecs,2),locs(resec_elecs,3),mark_size,...
    'rx','LineWidth',3);
scatter3(locs(:,1),locs(:,2),locs(:,3),mark_size,'k','LineWidth',3); 

view(-120,-11);


end