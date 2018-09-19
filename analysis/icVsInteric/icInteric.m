function icInteric(pt,whichPt,whichSz)

% output file name
[~,~,~,resultsFolder,~] = fileLocations;

% Get wij
xyChan = pt(whichPt).electrodeData.locs;
wij = getwij(xyChan,pt(whichPt).dmin);
nchs = length(pt(whichPt).channels);

% Get all sequences
seqs = pt(whichPt).sz(whichSz).seq_matrix;
seqs(seqs==0) = nan; % WHY ARE THERE ANY ZEROS?????

% Get seizure times
szTimes = [pt(whichPt).sz(whichSz).onset pt(whichPt).sz(whichSz).offset];

% Get the time of the first spike in each sequence
firstSpikes = min(seqs,[],1);

% separate sequences into ictal and interictal
seq_ic = seqs(:,firstSpikes >= szTimes(1) & firstSpikes <= szTimes(2));
seq_inter = seqs(:,~(firstSpikes >= szTimes(1) & firstSpikes <= szTimes(2)));


% get average recruitment latency for ictal and interictal sequences
[RL_ic,~] = getRL(seq_ic);
[RL_inter,~] = getRL(seq_inter);

% Get MI
MIstruct_ic = moranStats(RL_ic',wij,nchs);
MIstruct_inter = moranStats(RL_inter',wij,nchs); 


%% Test my new way of getting RL
%{
% confirm sequences are the same
seq_matrix = makeSeqMatrix(pt(whichPt).sz(whichSz).data.sequences,nchs);

% Old way
[recruitmentLatencySingle,spikeCount] = ...
    getRecruitmentLatency(pt(whichPt).sz(whichSz).data.sequences,xyChan);
[avgRecruitmentLat,I,MI] = getSpatialOrg(recruitmentLatencySingle,xyChan,pt(whichPt).dmin);
easyPlot(avgRecruitmentLat,xyChan);

% New way
[RL,~] = getRL(seqs);
easyPlot(RL,xyChan);
%}

end