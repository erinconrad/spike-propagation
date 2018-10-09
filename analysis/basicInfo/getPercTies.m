function pt = getPercTies(pt,whichPts)


for whichPt = whichPts
    for whichSz = 1:length(pt(whichPt).sz)
        if isfield(pt(whichPt).sz(whichSz),'data') == 0, continue, end
        fprintf('%s seizure %d has %1.2f percent ties\n',pt(whichPt).name,...
            whichSz,pt(whichPt).sz(whichSz).data.discarded.totalPercTies);
        
    
    end
end