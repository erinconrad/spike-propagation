function [all_seq_cat,all_times,seq_all,trackingNo] = divideIntoSzChunksGen(pt,whichPt)

clip_time = 600;

if isfield(pt(whichPt),'seq_matrix') == 0

% reinitialize the sequences
seq_all = {};

% just look at first seizure
sequences = pt(whichPt).sz(1).seq_matrix;
sequences(sequences==0) = nan; % WHY ARE THERE ANY ZEROS?????
firstSpikes = min(sequences,[],1);  

% remove sequences occuring during the seizure or ten minutes before or
% after the seizure
sz = [pt(whichPt).sz(1).onset - clip_time , pt(whichPt).sz(1).offset + clip_time];
sequences(:,firstSpikes>=sz(1) & firstSpikes<=sz(2)) = [];
seq_all{1} = sequences;


% Loop through the other seizures
for j = 2:length(pt(whichPt).sz)
    szTimeLast = pt(whichPt).sz(j-1).onset;
    szTime = pt(whichPt).sz(j).onset;
    totalTime = pt(whichPt).sz(j).runTimes(end,2) - pt(whichPt).sz(j).runTimes(1,1);
    
    sequences = pt(whichPt).sz(j).seq_matrix;
    sequences(sequences==0) = nan; % WHY ARE THERE ANY ZEROS?????

    
    if szTime - szTimeLast > totalTime
        % If the seizure time is more than 24 hours after the last seizure,
        % put this in a new chunk for plotting purposes
        seq_all{end+1} = sequences;
    else
        % if it is less than 24 hours after the last seizure, then need to
        % add any spikes that occured after the last seizure run time to
        % this new chunk (and ignore spikes from this seizure that would
        % have already been captured in the first seizure run)
        keepAfter = pt(whichPt).sz(j-1).runTimes(end,2);
        firstSpikes = min(sequences,[],1);  
        seq_to_keep = sequences(:,firstSpikes > keepAfter);

        % remove sequences occuring during the seizure
        sz = [pt(whichPt).sz(j).onset - clip_time pt(whichPt).sz(j).offset+ clip_time];
        firstSpikes = min(seq_to_keep,[],1);  
        seq_to_keep(:,firstSpikes>=sz(1) & firstSpikes<=sz(2)) = [];
        
        seq_all{end} = [seq_all{end},seq_to_keep];
    end

end



else
    
    
seq_all = cell(size(pt(whichPt).allTimes,1),1);
all_seq_cat = pt(whichPt).seq_matrix;
firstSpikes = min(all_seq_cat,[],1);

% Clip seizure times
szTimes = zeros(length(pt(whichPt).sz),2);
for j = 1:length(pt(whichPt).sz)
    szTimes(j,:) = [pt(whichPt).sz(j).onset pt(whichPt).sz(j).offset];
end

inSz = any(firstSpikes >= szTimes(:,1) & firstSpikes <= szTimes(:,2));

all_seq_cat(:,inSz) = [];
firstSpikes(inSz) = [];
    
% Split into chunks
for i = 1:size(pt(whichPt).allTimes,1)
    seq_all{i} = all_seq_cat(:,firstSpikes >= pt(whichPt).allTimes(i,1) & firstSpikes<= pt(whichPt).allTimes(i,2));
end

end

%% Recombine all sequences for the patient
all_seq_cat = [];
trackingNo = [];
for i = 1:length(seq_all)
    seq = seq_all{i};
    all_seq_cat = [all_seq_cat,seq];
    trackingNo = [trackingNo, i*ones(1,size(seq,2))];
end

all_times = min(all_seq_cat,[],1);




end