function aes_brain_colors(pt)

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
other_file_out = [outputFolder,'brain_colors'];


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
view(-120,-11);

end