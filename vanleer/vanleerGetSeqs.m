function seqs = vanleerGetSeqs(vanleer)

% The goal here is to make sequences that will be in the same format as my
% normal sequences

%% Parameters
minFrac = 0.1; %

times = vanleer.times;
rms = vanleer.rms;

%% Put in format of n_chs x n_spikes
rms = rms';
times = times';

seqs = nan(size(times));

% Loop through spikes
for i = 1:size(rms,2)
    
    % Find the indices of channels where the rms is above a required
    % minimum for that spike
    bigEnough = rms(:,i) > minFrac*max(rms(:,i));
    
    rms(~bigEnough,i) = nan;
    times(~bigEnough,i) = nan;
    
    seqs(:,i) = times(:,i);
    
end




end