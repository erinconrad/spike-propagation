clear

%% File names
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
%electrodeFile = [electrodeFolder,csvFile];
p1 = genpath(scriptFolder);
addpath(p1);
gdfFolder = 'gdf/';
chLocationsFolder = 'chLocations/';
ptWithSeq = 'ptWithSeq.mat';
finalPt = 'finalPt.mat';
ptWithIctal = 'ptWithIctal.mat';

load([resultsFolder,finalPt]);

% loop through patients
for i = 1:length(pt)

    % loop through seizures
    for j = 1:length(pt(i).sz)
       sztimes =  [pt(i).sz(j).onset,pt(i).sz(j).offset];
       % Loop through the blocks
       for k = 1:length(pt(i).sz(j).blockRL)
           blocktimes = pt(i).sz(j).blockRL(k).times + pt(i).sz(j).runTimes(1,1);
           
           pt(i).sz(j).blockRL(k).adjustedTimes = blocktimes;
           
           % if the start of the block is before the end of the seizure,
           % and the end of the block is after the start of the seizure,
           % then the block includes some seizure time
           if blocktimes(1) < sztimes(2) && blocktimes(2) > sztimes(1)
               pt(i).sz(j).blockRL(k).ictal = 1;
           else
               pt(i).sz(j).blockRL(k).ictal = 0;
           end
        
       end
    end

end


save([resultsFolder,ptWithIctal],'pt');