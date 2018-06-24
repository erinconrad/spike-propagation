function ictal = decideIfIctal(sequences,seizureTimes)
%{

This function takes a list of sequences and a list of seizure times and
returns an array of length equal to the number of sequences stating for
each sequence whether it is "ictal" or not. I call a sequence "ictal" if
the FIRST spike in the sequence occurs during a seizure.

%}

% number of sequences
nseq = size(sequences,2)/2;

% number of seizures
nseizures = size(seizureTimes,1);


% initialize output variables
ictal.ictal_idx = zeros(nseq,1);
ictal.seq_ictal = [];
ictal.seq_interictal = [];

% Loop through sequences
for s = 1:nseq
    
   % Get the time column and channel column for each sequence
   col = s*2;
   seqtime = sequences(:,col);
   seqch = sequences(:,col-1);
   
   % Loop through seizures
   for sz = 1:nseizures
       
       % Get seizure onset and offset times
       szTimes = seizureTimes(sz,:);
       
        % If the first spike in the sequence falls between the seizure
        % start and end time
        if seqtime(1) > szTimes(1) && seqtime(1) < szTimes(2)

           % it is an ictal sequence
           ictal.ictal_idx(s) = 1;
           
           % add the sequence times and channels to the list of ictal sequences
           ictal.seq_ictal = [ictal.seq_ictal,[seqch,seqtime]];
           
        end

   end
   
   % If, once you've looped through all the seizures, the ictal index is
   % still 0, then it's an interictal sequences
   if ictal.ictal_idx(s) == 0
      ictal.seq_interictal = [ictal.seq_interictal,[seqch,seqtime]];
       
   end
    
end




end



