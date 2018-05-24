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
finalPt = 'icvsinterict.mat';

load([resultsFolder,'ptStructs/',ptWithSeq]);

for i = 1:length(pt)
   for j = 1:length(pt(i).sz)
       
       if isfield(pt(i).sz(j),'data') == 0
          continue 
       end
       if isempty(pt(i).sz(j).data) == 1
           continue
       end
       
       szTimes = [pt(i).sz(j).onset-pt(i).sz(j).runTimes(1,1),...
           pt(i).sz(j).offset-pt(i).sz(j).runTimes(1,1)];
       
       nseq = pt(i).sz(j).stats.nseqs;
       chandata = pt(i).sz(j).data.xyChan;
       
       sequences = pt(i).sz(j).data.sequences;
       ictal = zeros(nseq,1);
       seq_ictal = [];
       seq_interictal = [];
       
       for s = 1:nseq
           col = s*2;
           seqtime = sequences(:,col);
           seqch = sequences(:,col-1);
           if seqtime(1) > szTimes(1) && seqtime(1) < szTimes(2)
               ictal(s) = 1;
               seq_ictal = [seq_ictal,[seqch,seqtime]];
           else
               seq_interictal = [seq_interictal,[seqch,seqtime]];
           end
           
       end
       
       pt(i).sz(j).data.ictal_idx = ictal;
       pt(i).sz(j).data.seq_ictal = seq_ictal;
       pt(i).sz(j).data.seq_interictal = seq_interictal;
       
       %% Get spatial org for the entire seizure period
 
       [rls_ictal,~] =getRecruitmentLatency(seq_ictal,chandata);
       
       [rl_ictal,so_ictal] = getSpatialOrg(rls_ictal,chandata,dmin);
       
       [rls_interictal,~] =getRecruitmentLatency(seq_interictal,chandata);
       
       [rl_interictal,so_interictal] = getSpatialOrg(rls_interictal,chandata,dmin);
         
       pt(i).sz(j).ic.rl_ictal = rl_ictal;
       pt(i).sz(j).ic.rl_interictal = rl_interictal;
       pt(i).sz(j).ic.so_ictal =so_ictal;
       pt(i).sz(j).ic.so_interictal = so_interictal;
       
       
       %% Now shuffle the ictal and interictal identities
       n_ic = length(find(ictal == 1));
       n_interic = length(find(ictal == 0));
       nperm = 1e3;
       all_so_diff = zeros(1,nperm);
       
       for iperm = 1:nperm
           if mod(iperm,10) == 0
               fprintf('Doing iteration %d\n',iperm)
               
           end
           
           all_seq = 1:nseq;
           new_ic = randperm(nseq,n_ic);
           new_interic = all_seq(~new_ic);
           
           for s = 1:nseq
               col = s*2;
               seqtime = sequences(:,col);
               seqch = sequences(:,col-1);
               if ismember(new_ic,s) == 1
                   seq_ictal = [seq_ictal,[seqch,seqtime]];
               else
                   seq_interictal = [seq_interictal,[seqch,seqtime]];
               end

           end
           
           [rls_ictal,~] =getRecruitmentLatency(seq_ictal,chandata);
       
           [rl_ictal,so_ictal] = getSpatialOrg(rls_ictal,chandata,dmin);

           [rls_interictal,~] =getRecruitmentLatency(seq_interictal,chandata);

           [rl_interictal,so_interictal] = getSpatialOrg(rls_interictal,chandata,dmin);
           
           all_so_diff(iperm) = so_interictal-so_ictal;
           
           
       end
       all_so_diff = sort(all_so_diff);
       all_so_diff(nperm*.05/2)
       pt(i).sz(j).ic.all_so_diff = all_so_diff;
       
   end
   
end

save([resultsFolder,'ptStructs/',finalPt],'pt');