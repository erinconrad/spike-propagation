function ptOld = mergeStructs(ptOld,ptTemp,whichPts)

for whichPt = whichPts

    ptOld(whichPt).thresh = ptTemp(whichPt).thresh;
    ptOld(whichPt).sz = ptTemp(whichPt).sz;
    ptOld(whichPt).allTimes = ptTemp(whichPt).allTimes;
    ptOld(whichPt).runTimes = ptTemp(whichPt).runTimes;
    ptOld(whichPt).chunkFiles = ptTemp(whichPt).chunkFiles;
    ptOld(whichPt).electrodeData = ptTemp(whichPt).electrodeData;
    ptOld(whichPt).channels = ptTemp(whichPt).channels;
    ptOld(whichPt).stats = ptTemp(whichPt).stats;
    ptOld(whichPt).data = ptTemp(whichPt).data;
    ptOld(whichPt).seq_matrix = ptTemp(whichPt).seq_matrix;

end

end