function info = getVectors(sequences,nseq,ictal,chanData)

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

   seqtime = sequences(1,col);

    % For each sequence, get the first and last spike in the sequence
    firstCh = realChs(1);
    lastCh = realChs(end);

    % Get the locations of the first and last spike, and calculate a direction
    % vector
    firstChLoc = chanData(firstCh,2:4);
    lastChLoc = chanData(lastCh,2:4);
    vector = lastChLoc -  firstChLoc;

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