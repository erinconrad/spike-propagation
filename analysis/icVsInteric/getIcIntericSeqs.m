function pt = getIcIntericSeqs(pt,whichPts)

for whichPt = whichPts
    
    fprintf('%s\n',pt(whichPt).name);
    
    for whichSz = 1:length(pt(whichPt).sz)
        
        if isfield(pt(whichPt).sz(whichSz),'seq_matrix') == 0
            continue
        end

        % Get all sequences
        seqs = pt(whichPt).sz(whichSz).seq_matrix;
        seqs(seqs==0) = nan; % WHY ARE THERE ANY ZEROS?????

        % Get seizure times
        szTimes = [pt(whichPt).sz(whichSz).onset pt(whichPt).sz(whichSz).offset];

        % Get the time of the first spike in each sequence
        firstSpikes = min(seqs,[],1);

        % separate sequences into ictal and interictal
        pt(whichPt).sz(whichSz).icinter.seq_ic = seqs(:,firstSpikes >= szTimes(1) & firstSpikes <= szTimes(2));
        pt(whichPt).sz(whichSz).icinter.seq_inter = seqs(:,~(firstSpikes >= szTimes(1) & firstSpikes <= szTimes(2)));
        pt(whichPt).sz(whichSz).icinter.icIdx = firstSpikes >= szTimes(1) & firstSpikes <= szTimes(2);
        
        fprintf('%d interictal sequences, %d ictal seqs\n',...
            size(pt(whichPt).sz(whichSz).icinter.seq_inter,2),...
            size(pt(whichPt).sz(whichSz).icinter.seq_ic,2));
        
        fprintf('Seizure duration is %1.1f s\n',...
            pt(whichPt).sz(whichSz).offset - pt(whichPt).sz(whichSz).onset);

    end

end

end