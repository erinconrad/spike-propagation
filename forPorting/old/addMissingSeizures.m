function pt = addMissingSeizures(pt)

%% Info on missing seizures
whichPt = 18; %HUP107
oldNumSz = size(pt(whichPt).newSzTimes,1);
missSz = [...
    721544.34 721569.68;...
    726981.46 726996.62;...
    730626.53 730640.62;...
    738893.83 738984.92;...
    746099.73 746116.84;...
    752517.47 752534.30];

%% Add seizure times

pt(whichPt).newSzTimes = [pt(whichPt).newSzTimes;missSz];

%% Add sz structures

for i = 1:size(missSz,1)
    pt(whichPt).sz(oldNumSz+i).onset = missSz(i,1);
    pt(whichPt).sz(oldNumSz+i).offset = missSz(i,2);
    pt(whichPt).sz(oldNumSz+i).electrodes = {};
    pt(whichPt).sz(oldNumSz+i).chs = [];
    
end

%% Expand all times
pt(whichPt).allTimes(2) = missSz(end,1) + 12*3600;

%% Expand run Times
chunkTime = 2000;
pt(whichPt).runTimes = [];
pt(whichPt).chunkFiles = {};
for j = 1:size(pt(whichPt).allTimes,1)
    initialTime = pt(whichPt).allTimes(j,1);
    finalTime = pt(whichPt).allTimes(j,2);
    totalTime = finalTime - initialTime;

    nchunks = ceil(totalTime/chunkTime);
    
    for k = 1:nchunks
       startTime = initialTime + (k-1)*chunkTime;
       endTime = min(finalTime, startTime + chunkTime);
       pt(whichPt).runTimes = [pt(whichPt).runTimes;startTime, endTime];

       % define the output file
       pt(whichPt).chunkFiles = [pt(whichPt).chunkFiles;...
           [pt(whichPt).name,'_subset_',sprintf('%d',j),'_times_',sprintf('%d',startTime),...
           '-',sprintf('%d',endTime),'.mat']];


   end
    
end

end