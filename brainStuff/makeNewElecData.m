function A = makeNewElecData(pt,whichPt)

elecFolder = '/Users/erinconrad/Desktop/residency stuff/R25/actual work/data/transformedElectrodes/';
%brainFolder = '/Users/erinconrad/Desktop/residency stuff/R25/actual work/data/brains/';

name = pt(whichPt).name;
names = dir([elecFolder,name,'*']);
fname = names.name;
fileID = fopen([elecFolder,fname]);
out=textscan(fileID, '%s', 'whitespace',',');
out = out{1};
nchs = length(out)/5;

for i = 1:nchs
    electrodeData.electrodes(i).x = str2double(out{(i-1)*5+1});
    electrodeData.electrodes(i).y = str2double(out{(i-1)*5+2});
    electrodeData.electrodes(i).z = str2double(out{(i-1)*5+3});
    electrodeData.electrodes(i).xyz = [electrodeData.electrodes(i).x,...
       electrodeData.electrodes(i).y,electrodeData.electrodes(i).z];
    electrodeData.electrodes(i).name = out{(i-1)*5+4};
end

%% Align this with the original electrode data
% Number of electrodes different because my main electrode data ignores
% electrodes that have no ieeg data
newData(whichPt).locs = [];
for i = 1:length(pt(whichPt).electrodeData.electrodes)
    for j = 1:length(electrodeData.electrodes)
        if strcmp(pt(whichPt).electrodeData.electrodes(i).name,...
                electrodeData.electrodes(j).name) == 1
            newData(whichPt).electrodes(i).name = electrodeData.electrodes(j).name;
            newData(whichPt).electrodes(i).xyz = electrodeData.electrodes(j).xyz;
            newData(whichPt).locs = [newData(whichPt).locs;newData(whichPt).electrodes(i).xyz];
        end

    end

end

if length(pt(whichPt).electrodeData.electrodes) ~= length(newData(whichPt).electrodes)
    error('Number of electrodes not aligned\n');
end

newLocs = newData(whichPt).locs;
oldLocs = pt(whichPt).electrodeData.locs(:,2:4);
A = newLocs*pinv(oldLocs);

%{
offset = [-1.0029,2.3087,28.9465];


plotLocs = A*oldLocs - offset;
plotLocsOrig = newLocs - offset;

%% Open gifti file
giftiFolder = [brainFolder,pt(whichPt).name,'/'];
names = dir([giftiFolder,'*pial.gii']);
fname2 = names.name;
g = gifti([giftiFolder,fname2]);

figure
p = plotGIFTI(g);
hold on
scatter3(plotLocsOrig(:,1),plotLocsOrig(:,2),...
    plotLocsOrig(:,3),80,'k','filled');
alpha(p,0.4)

figure
p = plotGIFTI(g);
hold on
scatter3(plotLocs(:,1),plotLocs(:,2),plotLocs(:,3),80,'k','filled');
alpha(p,0.4)

close(fileID);
%}
    



end