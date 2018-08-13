function getStatsOnAWindow(P,pt,sz,window)

sequences = P(pt).sz(sz).seq_matrix;
nseqs = size(sequences,2);
chLocs = P(pt).sz(sz).data.xyChan(:,2:4);
nchs = length(P(pt).channels);

% output file name
[~,~,~,resultsFolder,~] = fileLocations;

totalTime = max(sequences(:,end))-min(sequences(:,1));
if window > totalTime
    error('Window is too big\n'); 
end

nchunks = ceil(totalTime/window);


%% what to get stats on in that window
indices = cell(nchunks,1);
seqFreqTotal = zeros(nchunks,1); 
seqFreqChs = zeros(nchunks,nchs);
meantimes = zeros(nchunks,1);
angles = cell(nchunks,1);
mu = nan(nchunks,1);
ul = nan(nchunks,1);
ll = nan(nchunks,1);

freqChsInfo.save = [resultsFolder,'analysis/spikeFreqByCh/',P(pt).name...
    ,'_',sprintf('%d',window),'swindow.gif'];

for i = 1:nchunks
    times = [(i-1)*window+min(sequences(:,1)),i*window+min(sequences(:,1))];
    meantimes(i) = (times(1)+times(2))/2;
    
    %% get the appropriate sequences in this time
    % just look at first spike
    firstSpikes = min(sequences,[],1);    
    correctSequences = ...
        sequences(:,find(firstSpikes >= times(1) & firstSpikes <= times(2)));
    
    indices{i} = find(firstSpikes >= times(1) & firstSpikes <= times(2));
    
    seqFreqTotal(i) = size(correctSequences,2)/window;
    seqFreqChs(i,:) = nansum(correctSequences,2)/window;
    freqChsInfo.title{i} = sprintf('Sequence frequency by channel at time %d s',meantimes(i));
    
    angles{i} = getAngles(correctSequences,P(pt).electrodeData);
    if isempty(angles{i}) == 0
        [mu(i) ul(i) ll(i)] = circ_mean(angles{i}*pi/180);
    end
    
end

%% Get sequence bins
bins = zeros(nseqs,1);
for i = 1:nchunks
   if isempty(indices{i}) == 1
       continue
   end
   whichSeq = indices{i};
   bins(whichSeq(1):whichSeq(end)) = i;
    
end

allAngles = getAngles(sequences,P(pt).electrodeData);

%% Plots

scatter(1:length(mu),mu,'k','filled')
hold on
scatter(1:length(mu),ul,'k')
scatter(1:length(mu),ll,'k')

[pval table] = circ_wwtest(allAngles*pi/180,bins)

%{
% Total sequence frequency over time
figure
plot(meantimes/3600,seqFreqTotal)
xlabel('Time (hr)');
ylabel('Sequence frequency');
title(sprintf('Sequence frequency over time, window %d s',window));
set(gca,'FontSize',20)

brainMovieOfAnything(seqFreqChs,chLocs,freqChsInfo)
%}

end