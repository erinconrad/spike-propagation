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
icIdx = firstSpikes >= szTimes(1) & firstSpikes <= szTimes(2);

% get spike counts per ch (only those that show up in a sequence)
spikes_ch_all = sum(~isnan(seqs),2);
spikes_ch_ic = sum(~isnan(seq_ic),2);
spikes_ch_inter = sum(~isnan(seq_inter),2);

% get average recruitment latency for ictal and interictal sequences
[RL,~] = getRL(seqs);
[RL_ic,~] = getRL(seq_ic);
[RL_inter,~] = getRL(seq_inter);

% Get MI
MIstruct_ic = moranStats(RL_ic',wij,nchs);
MIstruct_inter = moranStats(RL_inter',wij,nchs); 

% Get vectors
vec_all = getVectors2(seqs,xyChan);
vec_ic = getVectors2(seq_ic,xyChan);
vec_inter = getVectors2(seq_inter,xyChan);

% MANOVA comparing vectors for ictal and interictal sequences
[d,p,stats] = manova1(vec_all,icIdx,0.01);

% the next step will be to figure out how to combine it for multiple
% seizures and multiple patients with manova
% https://www.mathworks.com/help/stats/repeatedmeasuresmodel.manova.html


% Test that the vectors form more or less a normal distribution
%{
mean_vec = mean(vec_all,1);
std_vec = std(vec_all,1);
n_vecs = size(vec_all,1);
z = repmat(mean_vec,n_vecs,1) + randn(n_vecs,3).*std_vec;
scatter3(z(:,1),z(:,2),z(:,3),'g');
hold on
scatter3(vec_all(:,1),vec_all(:,2),vec_all(:,3),'b');
%}

% Compare the ictal and interictal vectors spatially
%{
figure
scatter3(vec_ic(:,1),vec_ic(:,2),vec_ic(:,3),'b')
hold on
scatter3(vec_inter(:,1),vec_inter(:,2),vec_inter(:,3),'r')

%}


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