function aes_build_network(pt)

% This is for HUP078
whichPt = 8;
offset = [-3 27 7.7653];
mark_size = 600;

white_elecs = [70 76 85];

%% Get stuff
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile,other] = fileLocations;
addpath(other.gifti)
addpath(genpath(scriptFolder))


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
scatter3(locs(white_elecs,1),locs(white_elecs,2),locs(white_elecs,3),...
    mark_size,'markerfacecolor',[1 1 1]);
scatter3(locs(:,1),locs(:,2),locs(:,3),mark_size,'k','LineWidth',3); 

for i = 2:length(white_elecs)
    plot3([locs(white_elecs(i-1),1) locs(white_elecs(i),1)],...
        [locs(white_elecs(i-1),2) locs(white_elecs(i),2)],...
        [locs(white_elecs(i-1),3) locs(white_elecs(i),3)],...
        'r','linewidth',3)
end
alpha(p, 0.3)


view(-120,-11);

end