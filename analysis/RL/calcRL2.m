% This script calculates the spatial organization for each block and for
% the whole seizure period.


clear

%% Parameters


%% File names
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
%electrodeFile = [electrodeFolder,csvFile];
p1 = genpath(scriptFolder);
addpath(p1);
ptWithSeq = 'ptWithSeq.mat';
finalPt = 'finalPt.mat';

load([resultsFolder,'ptStructs/',ptWithSeq]);


for i = 1:length(pt)
   for j = 1:length(pt(i).sz)
       
       if isfield(pt(i).sz(j),'data') == 0
          continue 
       end
       if isempty(pt(i).sz(j).data) == 1
           continue
       end
       
       nseq = pt(i).sz(j).stats.nseqs;
       
       % Calculate recruitment latencies and spatial orgs for each block
       sequences = pt(i).sz(j).data.sequences;
       
       [RLS,spikeCount] = ...
              getRecruitmentLatency(sequences,pt(i).sz(j).data.xyChan);
       
   end
end

fprintf('look\n');
