%{
This script separates the sequences into ictal sequences (those that occur
during a seizure) and interictal sequences (those that don't occur during a
seizure). It gets the average recruitment latencies for the ictal sequences
and the interictal sequences. It then sorts the latencies and calculates a
direction vector going from the early channels (early 50%) to the late
channels (late 50%). It then measures the angle between the direction
vectors for the interictal sequences and the ictal sequences.

If there were no difference between the ictal and the interictal data, then
the angle between them would be 0.


%}



clear

%% Bonus parameter
dmin = 15;
nperm = 1e3;
minSeq = 10;

%% File names
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
%electrodeFile = [electrodeFolder,csvFile];
p1 = genpath(scriptFolder);
addpath(p1);
ptWithSeq = 'ptWithSeq.mat';
finalPt = 'finalPt.mat';
newPt = 'icIntericAngle.mat';

load([resultsFolder,'ptStructs/',ptWithSeq]);

% Loop through patients
for i = 1%1:length(pt)
    
    % Loop through seizures
   for j = 1:length(pt(i).sz)
       
       if isfield(pt(i).sz(j),'data') == 0
          continue 
       end
       if isempty(pt(i).sz(j).data) == 1
           continue
       end
       
       % Get the seizure times
       szTimes = [pt(i).sz(j).onset,...
           pt(i).sz(j).offset];
       
       nseq = pt(i).sz(j).stats.nseqs;
       chandata = pt(i).sz(j).data.xyChan;
       
       % Get the sequences
       sequences = pt(i).sz(j).data.sequences;
       
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
       
       % The ictal and interictal sequences
       pt(i).sz(j).data.seq_ictal = seq_ictal;
       pt(i).sz(j).data.seq_interictal = seq_interictal;
       
       %% Get spatial org for ictal and interictal sequences
 
       [rls_ictal,~] =getRecruitmentLatency(seq_ictal,chandata);
       
       [rl_ictal,so_ictal] = getSpatialOrg(rls_ictal,chandata,dmin);
       
       [rls_interictal,~] =getRecruitmentLatency(seq_interictal,chandata);
       
       [rl_interictal,so_interictal] = getSpatialOrg(rls_interictal,chandata,dmin);
         
       pt(i).sz(j).ic.rl_ictal = rl_ictal;
       pt(i).sz(j).ic.rl_interictal = rl_interictal;
       pt(i).sz(j).ic.so_ictal =so_ictal;
       pt(i).sz(j).ic.so_interictal = so_interictal;
       pt(i).sz(j).ic.so_diff = so_interictal-so_ictal;
       
       %% Now getting angles
       
       % Sort the recruitment latencies for the ictal sequences
       [~,A] = sort(rl_ictal);
       
       % The first half are the earliest channels recruited. Take the mean
       % of their position
       earlychans = mean(chandata(A(1:length(A)/2),2:4));
       
       % The second half are the latest channels recruited. Take the mean
       % of their position
       latechans = mean(chandata(A(length(A)/2+1:length(A)),2:4));
       
       % The direction vector going from the early channels to the late
       % channles
       direction_ic = latechans-earlychans;
       
       pt(i).sz(j).ic.earlychans_ic = earlychans;
       pt(i).sz(j).ic.latechans_ic = latechans;
       pt(i).sz(j).ic.direction_ic = direction_ic;
       
       % Do the same for interictal data
       [~,A] = sort(rl_interictal);
       earlychans = mean(chandata(A(1:length(A)/2),2:4));
       latechans = mean(chandata(A(length(A)/2+1:length(A)),2:4));
       direction_interic = latechans-earlychans;
       
       % Get the angle in degrees between the ictal direction and the
       % interictal direction
       angle = acosd(dot(direction_ic,direction_interic)/norm(direction_ic)/norm(direction_interic));
       
       pt(i).sz(j).ic.angle = angle;
       pt(i).sz(j).ic.earlychans_interic = earlychans;
       pt(i).sz(j).ic.latechans_interic = latechans;
       pt(i).sz(j).ic.direction_interic = direction_interic;
       
       %% Now shuffle the ictal and interictal identities
       
       % Get the number of ictal and interictal sequences
       n_ic = length(find(ictal == 1));
       n_interic = length(find(ictal == 0));
       
       % Initialize the array for the difference measure (the main measure
       % for the permutation testing)
       angles_perm = zeros(1,nperm);
       
       % Loop through the permutations
       for iperm = 1:nperm
           if mod(iperm,100) == 0
               fprintf('Doing iteration %d\n',iperm)
           end
           
           % Get a list from 1 to nseq, indices of all the sequences
           all_seq = 1:nseq;
           
           % Randomly select from this list a number of "ictal" sequences,
           % the same number as the true number of ictal sequences
           new_ic = randperm(nseq,n_ic);
           
           % Set the logical indices of these "ictal" sequences to 1
           new_ic_logic = zeros(size(all_seq));
           new_ic_logic(new_ic) = 1;
           
           % Get the indices of the "interictal" sequences
           new_interic = all_seq(~new_ic_logic);
           
           % Initialize arrays for these sequences. They are as long
           % vertically as the sequences were for the original sequences,
           % and as wide as the number of ictal sequences x 2 and
           % interictal sequences x 2, respectively (one time column and
           % one channel column)
           seq_ictal = zeros(size(sequences,1),2*n_ic);
           seq_interictal = zeros(size(sequences,1),2*n_interic);
           
           % Loop through the sequences
           for s = 1:nseq
               
               % Get the channel and time columns
               col = s*2;
               seqtime = sequences(:,col);
               seqch = sequences(:,col-1);
               
               % if this sequence is one of the new "ictal" sequences
               if ismember(s,new_ic) == 1
                   
                   % Fill up the ictal array
                   seq_ictal(:,col-1:col) = [seqch,seqtime];
               else
                   
                   % If not, fill up the interictal array
                   seq_interictal(:,col-1:col) = [seqch,seqtime];
               end

           end
           
           
           
           [rls_ictal,~] =getRecruitmentLatency(seq_ictal,chandata);
       
           [rl_ictal,so_ictal] = getSpatialOrg(rls_ictal,chandata,dmin);

           [rls_interictal,~] =getRecruitmentLatency(seq_interictal,chandata);

           [rl_interictal,so_interictal] = getSpatialOrg(rls_interictal,chandata,dmin);
           
           [~,A] = sort(rl_ictal);
           earlychans = mean(chandata(A(1:length(A)/2),2:4));
           latechans = mean(chandata(A(length(A)/2+1:length(A)),2:4));
           direction_ic = latechans-earlychans;

           [~,A] = sort(rl_interictal);
           earlychans = mean(chandata(A(1:length(A)/2),2:4));
           latechans = mean(chandata(A(length(A)/2+1:length(A)),2:4));
           direction_interic = latechans-earlychans;
           angle = acosd(dot(direction_ic,direction_interic)/norm(direction_ic)/norm(direction_interic));
           
           angles_perm(iperm) = angle;
           
           
       end
       
       % These are the angles you get from the iterations of the
       % permutation test, sorted from smallest to biggest
       angles_perm = sort(angles_perm);
       pt(i).sz(j).ic.angles_perm = angles_perm;
       
       %% Get a p value
       
       % take the difference between the array of permutation test angles
       % and the true angle
       diffDiff = abs(pt(i).sz(j).ic.angle-angles_perm);
       
       % Find where the difference is smallest. This is essentially how
       % significant the angle is (if it's really big or small it will be
       % close to the top or bottom of the permutation test angles)
       [~,I] = min(diffDiff);
       
       % Need to multiply by 2 to get the actual p value
       p = 2*min(I/nperm,(1-I/nperm));
       pt(i).sz(j).ic.p = p;
       
   end
   
end

save([resultsFolder,'ptStructs/',newPt],'pt');