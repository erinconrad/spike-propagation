function numSeqs(pt,whichPts)

for whichPt = whichPts
    
    fprintf('%s\n',pt(whichPt).name);
        
    for i = 1:length(pt(whichPt).sz)
        fprintf('%d interictal sequences, %d ictal seqs\n',...
            size(pt(whichPt).sz(i).icinter.seq_inter,2),...
            size(pt(whichPt).sz(i).icinter.seq_ic,2));
        
        fprintf('Seizure duration is %1.1f s\n',...
            pt(whichPt).sz(i).offset - pt(whichPt).sz(i).onset);
        
    end


end

end