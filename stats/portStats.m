clear

% Get spike frequency, sequence frequency,
% interictal sequence frequency, percentage tied

%% File names
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
%electrodeFile = [electrodeFolder,csvFile];
p1 = genpath(scriptFolder);
addpath(p1);
ptWithIctal = 'ptWithIctal.mat';
ptWithStats = 'ptWithStats.mat';

allSpF = [];
allSeF = [];

load([resultsFolder,ptWithIctal]);


for i = 1:length(pt)
   for j = 1:length(pt(i).sz)
       if isfield(pt(i).sz(j),'stats') == 0
           continue
       end
       pt(i).sz(j).stats.fspikespers = pt(i).sz(j).stats.nspikes/...
           (pt(1).sz(1).runTimes(end,2)-pt(1).sz(1).runTimes(1,1));
       pt(i).sz(j).stats.fseqpermin = pt(i).sz(j).stats.nseqs/...
            (pt(1).sz(1).runTimes(end,2)-pt(1).sz(1).runTimes(1,1))*60;
       pt(i).sz(j).stats.percTies = pt(i).sz(j).data.discarded.totalPercTies;
       
       allSpF = [allSpF,pt(i).sz(j).stats.fspikespers];
       allSeF = [allSeF,pt(i).sz(j).stats.fseqpermin];
       
   end
end

save([resultsFolder,ptWithStats],'pt');