clear
dmin = 15;

%% File names
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
%electrodeFile = [electrodeFolder,csvFile];
p1 = genpath(scriptFolder);
addpath(p1);
finalPt = 'finalPt.mat';

load([resultsFolder,'ptStructs/',finalPt]);


for i = 1:length(pt)
   for j = 1:length(pt(i).sz)
       if isfield(pt(i).sz(j),'data') == 0
          continue 
       end
       if isempty(pt(i).sz(j).data) == 1
           continue
       end
       
       sequences = pt(i).sz(j).data.sequences;
       chandata = pt(i).sz(j).data.xyChan;
       nseq = pt(i).sz(j).stats.nseqs;
       
       szTimes = [pt(i).sz(j).onset-pt(i).sz(j).runTimes(1,1),...
           pt(i).sz(j).offset-pt(i).sz(j).runTimes(1,1)];
       
       %% ictal vs interictal
       
       ictal = zeros(nseq,1);
       seq_ictal = [];
       seq_interictal = [];
       
       % Loop through the sequences
       for s = 1:nseq
           
           % Get the time column and channel column for each sequence
           col = s*2;
           seqtime = sequences(:,col);
           seqch = sequences(:,col-1);
           
           % If the first spike in the sequence falls between the seizure
           % start and end time
           if seqtime(1) > szTimes(1) && seqtime(1) < szTimes(2)
               
               % it is an ictal sequence
               ictal(s) = 1;
               
               % add the sequence times and channels to the list of ictal sequences
               seq_ictal = [seq_ictal,[seqch,seqtime]];
           else
               
               % it's not an ictal sequence, add it to the list of
               % interictal sequences
               seq_interictal = [seq_interictal,[seqch,seqtime]];
           end
           
       end
       
       % the indices of ictal sequences
       pt(i).sz(j).data.ictal_idx = ictal;
       pt(i).sz(j).data.seq_ictal = seq_ictal;
       pt(i).sz(j).data.seq_interictal = seq_interictal;
       
       %% Get spatial org for the entire seizure period
 
       [rls_ictal,~] =getRecruitmentLatency(seq_ictal,chandata);
       
       [rl_ictal,so_ictal,MI_ictal] = getSpatialOrg(rls_ictal,chandata,dmin);
       
       [rls_interictal,~] =getRecruitmentLatency(seq_interictal,chandata);
       
       [rl_interictal,so_interictal,MI_interictal] = getSpatialOrg(rls_interictal,chandata,dmin);
         
       Z = (MI_ictal.I-MI_interictal.I)/sqrt(MI_ictal.V+MI_interictal.V);
       p = normcdf(-Z);
       pt(i).sz(j).ic.rl_ictal = rl_ictal;
       pt(i).sz(j).ic.rl_interictal = rl_interictal;
       
       [~,A] = sort(rl_ictal);
       earlychans = mean(chandata(A(1:length(A)/2),2:4));
       latechans = mean(chandata(A(length(A)/2+1:length(A)),2:4));
       direction_ic = latechans-earlychans;
       
       [~,A] = sort(rl_interictal);
       earlychans = mean(chandata(A(1:length(A)/2),2:4));
       latechans = mean(chandata(A(length(A)/2+1:length(A)),2:4));
       direction_interic = latechans-earlychans;
       angle = acosd(dot(direction_ic,direction_interic)/norm(direction_ic)/norm(direction_interic));
       
   end
   
end