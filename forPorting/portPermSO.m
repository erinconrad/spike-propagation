clear

%% Bonus parameter
dmin = 15;
nperm = 1e3;
minSeq = 10;
alpha = 0.05;

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
          fprintf('Doing permutation test for block %d of %d\n',t,length(pt(i).sz(j).blockRL));
           
          sIdx = pt(i).sz(j).blockRL(t).sIdx;
          newseq = [];
          % Get the columns from the sIdx
          for k = 1:length(sIdx)
              newseq = [newseq,sequences(:,k*2-1:k*2)];
          end
          
          %% Do it once regularly
          
          [recruitmentLatencySingle,spikeCount] = ...
              getRecruitmentLatency(newseq,pt(i).sz(j).data.xyChan);
          
          [pt(i).sz(j).blockRL(t).avgRecruitmentLat,...
              pt(i).sz(j).blockRL(t).spatialOrg] = ...
              getSpatialOrg(recruitmentLatencySingle,pt(i).sz(j).data.xyChan,dmin);
          
          spatialOrg = pt(i).sz(j).blockRL(t).spatialOrg;
          MIstat = zeros(nperm,1);
          chandata = pt(i).sz(j).data.xyChan;
          %% Loop through permutations
          for ip = 1:nperm
              
           
              
              perm = randperm(size(chandata,1));
              newchan = chandata(perm,:);
              
              [rl,~] = getRecruitmentLatency(newseq,newchan);
              [~,MIstat(ip)] = getSpatialOrg(rl,newchan,dmin);
          end
          
          MIstat = sort(MIstat);
          CI = [MIstat(max(round(alpha*nperm/2),1)),...
              MIstat(min(round(length(MIstat)-alpha*nperm/2),length(MIstat)))];
          soDiff = abs(spatialOrg-MIstat);
          [~,I] = min(soDiff);
          p = 2*min(I/nperm,(1-I/nperm));
          pt(i).sz(j).blockRL(t).CI95 = CI;
          pt(i).sz(j).blockRL(t).p = p;
          
          % re-define the spatial org as nan if not enough sequences
          %{
          if length(sIdx) < minSeq
              pt(i).sz(j).blockRL(t).oldSpatialOrg = pt(i).sz(j).blockRL(t).spatialOrg;
              pt(i).sz(j).blockRL(t).spatialOrg = nan;
              pt(i).sz(j).blockRL(t).oldCI = pt(i).sz(j).blockRL(t).CI95;
              pt(i).sz(j).blockRL(t).CI95 = [nan nan];
          end
          %}
           
       end
       
       %% Get spatial org for the entire seizure period
       
       % Do it once regularly
       [recruitmentLatencySingleAll,spikeCountAll] =...
           getRecruitmentLatency(sequences,chandata);
       
       [pt(i).sz(j).avgRecruitmentLatAll,...
           pt(i).sz(j).spatialOrgAll] = ...
           getSpatialOrg(recruitmentLatencySingleAll,...
             chandata,dmin);
         
        % Permutation test
        MIstatAll = zeros(nperm,1);
        for ip = 1:nperm
          perm = randperm(size(chandata,1));
          newchan = chandata(perm,:);

          [rl,~] = getRecruitmentLatency(sequences,newchan);
          [~,MIstatAll(ip)] = getSpatialOrg(rl,newchan,dmin);
        end
        
        MIstatAll = sort(MIstatAll);
        CIAll = [MIstatAll(max(1,round(alpha*nperm/2))),MIstatAll(min(length(MIstat),round(length(MIstat)-alpha*nperm/2)))];
        soDiffAll = abs(spatialOrg-MIstatAll);
        [~,I] = min(soDiffAll);
        pAll = 2*min(I/nperm,(1-I/nperm));
        pt(i).sz(j).CI95All = CIAll;
        pt(i).sz(j).pAll = pAll;
   
   end
 
end

save([resultsFolder,'ptStructs/',finalPt],'pt');