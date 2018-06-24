function info = getVectors(sequences,nseq,ictal,chanData,vectorMethod)

nictal = sum(ictal);
ninterictal = length(ictal) - nictal;

ic_vectors = zeros(nictal,3);
interic_vectors = zeros(ninterictal,3);
ic_lineLength = zeros(nictal,1);
interic_linelength = zeros(ninterictal,1);

ic_idx = 1;
interic_idx = 1;

for s = 1:nseq
   col = s*2;
   seqchs = sequences(:,col-1);
   realChs = seqchs(seqchs~=0);

   seqtimes = sequences(:,col);

    % For each sequence, get the first and last spike in the sequence
    firstCh = realChs(1);
    lastCh = realChs(end);

    if vectorMethod == 1
        % Get the locations of the first and last spike, and calculate a direction
        % vector
        firstChLoc = chanData(firstCh,2:4);
        lastChLoc = chanData(lastCh,2:4);
        vector = lastChLoc -  firstChLoc;
    elseif vectorMethod == 2
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
        
        vector = secondHalfChLocsMean - firstHalfChLocsMean;
        
    end

    % Get a measure of line length (how much it jumps around)
    lineLength = 0;
    for k = 1:length(realChs)-1
        lineLength = lineLength + norm(chanData(realChs(k+1),2:4) - ...
            chanData(realChs(k),2:4)); 
    end
    lineLength = lineLength/length(realChs);

    % is it ictal or interictal
    if ictal(s) == 1
        ic_vectors(ic_idx,:) = vector;
        ic_lineLength(ic_idx) = lineLength;
        ic_idx = ic_idx + 1;
    else
        interic_vectors(interic_idx,:) = vector;
        interic_lineLength(interic_idx) = lineLength;
        interic_idx = interic_idx + 1;
    end


end

info.vectors.interic_vectors = interic_vectors;
info.vectors.ic_vectors = ic_vectors;
info.lineLength.ic_lineLength = ic_lineLength;
info.lineLength.interic_lineLength =interic_lineLength;

end