function aes_source_loc(pt)
% This is for HUP078
whichPt = 8;
mark_size = 600;

%% Get stuff
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile,other] = fileLocations;
addpath(other.gifti)
addpath(genpath(scriptFolder))
outputFolder = [resultsFolder,'pretty_plots/aes_gif/'];
mkdir(outputFolder)
other_file_out = [outputFolder,'brain_colors'];


locs = pt(whichPt).electrodeData.locs(:,2:4);

fig = figure;
set(fig,'position',[10 10 1000 800])
set(gcf,'color','white');

x = locs(:,1);
y = locs(:,2);
z = locs(:,3);

% Define a max at a certain point in the figure
center_point = [162,140,113];
max_dist = 30;
col_lin = zeros(size(x,1),1);
for i = 1:length(col_lin)
    dist = vecnorm(locs(i,:)-center_point);
    if dist >= max_dist
        col_lin(i) = 0;
    elseif dist <1
        col_lin(i) = 1;
    else
        col_lin(i) = 1./dist;
    end
end

%col_lin = arrayfun(@(x) max(log(x),-10),col_lin);
col_lin = arrayfun(@(x) (log((x))),col_lin);
%col_lin = sqrt(col_lin);

scatter3(x,y,z,mark_size,col_lin,'filled');
view(173,-45)
grid off
xticklabels([])
yticklabels([])
zticklabels([])
axis off


end