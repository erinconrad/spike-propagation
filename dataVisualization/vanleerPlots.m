clear

%% File names
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
%electrodeFile = [electrodeFolder,csvFile];
p1 = genpath(scriptFolder);
addpath(p1);
ptWithFs = 'finalPt.mat';
gdfFolder = [resultsFolder,'gdf/'];
chLocationsFolder = 'chLocations/';
ptVanleer = 'ptVanleer.mat';

load([resultsFolder,'ptStructs/',ptVanleer]);


gold_delay = pt(1).vanleer.gold_delay;
gold_rms = pt(1).vanleer.gold_rms;

for t = 1:10
   plotDelays(gold_delay(t,:),gold_rms(t,:),pt(1).electrodeData.locs); 

end