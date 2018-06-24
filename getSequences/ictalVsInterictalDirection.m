function ictalVsInterictalDirection(pt,whichPt,whichSz)

% which method to use for calculating vectors
method = 2;


% Which sequences
sequences = pt(whichPt).sz(whichSz).data.sequences;

%% Get ALL seizure times
nsz = length(pt(whichPt).sz);
seizureTimes = zeros(nsz,2);

for sz = 1:nsz
    seizureTimes(sz,:) = [pt(whichPt).sz(sz).onset,pt(whichPt).sz(sz).offset];
end

%% Get ictal versus interictal
ictal = decideIfIctal(sequences,seizureTimes);


%% Get directions
vectors  = getVectorsGeneral(sequences,method);

%% Get average ictal and interictal vectors
vector_ic = mean(vectors(ictal.ictal_idx));

end