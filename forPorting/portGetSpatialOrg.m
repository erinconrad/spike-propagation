% This script calculates the spatial organization for each block and for
% the whole seizure period.


clear

%% Bonus parameter
dmin = 15;
nperm = 1e4;
minSeq = 10;

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
       
       times = pt(i).sz(j).timesRL;
       nseq = pt(i).sz(j).stats.nseqs;
       
       for t = 1:size(times,1)
          pt(i).sz(j).blockRL(t).sIdx = [];
          
          % find what sequences are in the block
          for s = 1:nseq
              
              col = s*2;
              seqtime = pt(i).sz(j).data.sequences(1,col);
              if seqtime >= times(t,1) && seqtime <= times(t,2)
                  pt(i).sz(j).blockRL(t).sIdx = ...
                      [pt(i).sz(j).blockRL(t).sIdx,s];
              elseif seqtime < times(t,1)
                  continue
              elseif seqtime > times(t,2)
                  break
              end
          end
           
       end
       
       % Calculate recruitment latencies and spatial orgs for each block
       sequences = pt(i).sz(j).data.sequences;
       for t = 1:length(pt(i).sz(j).blockRL)
          sIdx = pt(i).sz(j).blockRL(t).sIdx;
          newseq = [];
          % Get the columns from the sIdx
          for k = 1:length(sIdx)
              newseq = [newseq,sequences(:,k*2-1:k*2)];
          end
          
          [recruitmentLatencySingle,spikeCount] = ...
              getRecruitmentLatency(newseq,pt(i).sz(j).data.xyChan);
          
          [pt(i).sz(j).blockRL(t).avgRecruitmentLat,...
              pt(i).sz(j).blockRL(t).spatialOrg,~] = ...
              getSpatialOrg(recruitmentLatencySingle,pt(i).sz(j).data.xyChan,dmin);
          
          % re-define the spatial org as nan if not enough sequences
          if length(sIdx) < minSeq
              pt(i).sz(j).blockRL(t).oldSpatialOrg = pt(i).sz(j).blockRL(t).spatialOrg;
              pt(i).sz(j).blockRL(t).spatialOrg = nan;
          end
           
       end
       
       % Get spatial org for the entire seizure period
       [recruitmentLatencySingleAll,spikeCountAll] =...
           getRecruitmentLatency(sequences,pt(i).sz(j).data.xyChan);
            
       [pt(i).sz(j).avgRecruitmentLatAll,...
           pt(i).sz(j).spatialOrgAll,pt(i).sz(j).MI] = ...
           getSpatialOrg(recruitmentLatencySingleAll,...
             pt(i).sz(j).data.xyChan,dmin);
   
   end
 
end

save([resultsFolder,'ptStructs/',finalPt],'pt');