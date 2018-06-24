function vectors = getVectorsGeneral(sequences,method)

%{

This function takes a list of sequences and returns direction vectors
summarizing the direction of the spike sequence

%}

% number of sequences
nseq = size(sequences,2)/2;

% initialize vectors
vectors = zeros(nseq,3);

% Loop through sequences
for s = 1:nseq
    col = s*2;
    seqchs = sequences(:,col-1);
    realChs = seqchs(seqchs~=0);
    seqtimes = sequences(:,col);
   
   % For each sequence, get the first and last spike in the sequence
    firstCh = realChs(1);
    lastCh = realChs(end);

    if method == 1
        % Get the locations of the first and last spike, and calculate a direction
        % vector
        firstChLoc = chanData(firstCh,2:4);
        lastChLoc = chanData(lastCh,2:4);
        vectors(s,:) = lastChLoc -  firstChLoc;
    elseif method == 2
        % Get the location of the mean of locations of the first half of
        % spikes in the sequence, and get the location of the mean of
        % locations of the 2nd half of spikes in the sequence. Calculate a
        % direction vector between those two means.
        
        % Get identities of channels
        firstHalfChs = realChs(1:floor(length(realChs)/2));
        secondHalfChs = realChs(ceil(length(realChs)/2):end);
        
        % Get mean locations of these channels
        firstHalfChLocsMean = mean(chanData(firstHalfChs,2:4));
        secondHalfChLocsMean = mean(chanData(secondHalfChs,2:4));
        
        vectors(s,:) = secondHalfChLocsMean - firstHalfChLocsMean;
        
    end
    
end

end