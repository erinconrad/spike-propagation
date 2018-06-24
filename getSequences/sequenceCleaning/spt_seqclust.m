%{

    Cluster spike sequences using clue-aware trajectory comparison

    Author: Sam Tomlinson, Aug 2016

%}

function seq_track = spt_seqclust(xyChan,allseqs,ss_thresh,tt_thresh)


% Generate files for clueGraph algorithm
[xlocs,ylocs,zlocs,chanlocs,tstamps] = spt_filegen(allseqs,xyChan);

% Cluster sequences using clueGraph method
[clueGraph] = spt_clueaware(ss_thresh,tt_thresh,xlocs,ylocs,zlocs,chanlocs,tstamps);

% Clean sequences
[clueGraph_clean,seq_track] = spt_cleanseqs(clueGraph);


end
    


















