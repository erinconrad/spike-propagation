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

% the indices of the sequences that fall into the chunk
indices = cell(nchunks,1);

% the sequence frequency in that chunk
seqFreqTotal = zeros(nchunks,1); 

% sequence frequency in each channel in that chunk
seqFreqChs = zeros(nchs,nchunks);

% mean time of the chunk
meantimes = zeros(nchunks,1);

% angles of all of the sequences that fall in that chunk
angles = cell(nchunks,1);

% mean angle, upper 95% CI, and lower 95% CI of the angles of the sequences
% in theat chunk
mu = nan(nchunks,1);
ul = nan(nchunks,1);
ll = nan(nchunks,1);

freqChsInfo.save = [resultsFolder,'analysis/spikeFreqByCh/',P(pt).name...
    ,'_',sprintf('%d',window),'swindow.gif'];

% Loop through the chunks
for i = 1:nchunks
    
    % get the correct times
    times = [(i-1)*window+min(sequences(:,1)),i*window+min(sequences(:,1))];
    meantimes(i) = (times(1)+times(2))/2;
    
    %% get the appropriate sequences in this time
    
    % just look at first spike
    firstSpikes = min(sequences,[],1);   
    
    % the correct sequences are those 
    correctSequences = ...
        sequences(:,find(firstSpikes >= times(1) & firstSpikes <= times(2)));
    
    indices{i} = find(firstSpikes >= times(1) & firstSpikes <= times(2));
    
    % Get the number of sequences in the window divided by the duration of
    % the window in seconds to get the sequence frequency in that window
    seqFreqTotal(i) = size(correctSequences,2)/window;
    
    % Count up all the non nans for each channel across the correct
    % sequences to get the number of sequences that channel is involved in
    % in the window
    seqFreqChs(:,i) = sum(~isnan(correctSequences),2)/window;
    freqChsInfo.title{i} = sprintf('Sequence frequency by channel at time %d s',meantimes(i));
    
    % need to check this
    %{
    angles{i} = getAngles(correctSequences,P(pt).electrodeData);
    if isempty(angles{i}) == 0
        [mu(i) ul(i) ll(i)] = circ_mean(angles{i}*pi/180);
    end
    %}
    
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

%allAngles = getAngles(sequences,P(pt).electrodeData);

%% Plots

%{
brainImageOfAnything(mean(seqFreqChs,2),chLocs,0);

scatter(1:length(mu),mu,'k','filled')
hold on
scatter(1:length(mu),ul,'k')
scatter(1:length(mu),ll,'k')

[pval table] = circ_wwtest(allAngles*pi/180,bins)
%}

brainMovieOfAnything(seqFreqChs',chLocs,freqChsInfo)

%{
% Total sequence frequency over time
figure
plot(meantimes/3600,seqFreqTotal)
xlabel('Time (hr)');
ylabel('Sequence frequency');
title(sprintf('Sequence frequency over time, window %d s',window));
set(gca,'FontSize',20)


%}

end