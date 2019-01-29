function newTimes = makeNonIntersectingTimeRanges(originalTimes,excludedTimes)

% How much give for something to be considered intersecting?
foregiveness = 1; %1 second


newTimes = [];

startTime = originalTimes(1,1);
whichOrig = 1;
whichExc = 1;

while 1
    
    currRange = [startTime originalTimes(whichOrig,2)];
    % flag if time period i doesn't intersect with any excluded time
    
    intersect = doTimeRangesIntersectForgiving(...
        currRange,excludedTimes(whichExc,:),0);
    if intersect == 1
        newTimes = [newTimes;startTime excludedTimes(whichExc,1)];
        startTime = excludedTimes(whichExc,2);
        whichExc = whichExc + 1;
    else
        whichExc = whichExc + 1;
    end
        
    if whichExc == size(excludedTimes,1) + 1
        
        % add the remaining times
        newTimes = [newTimes;startTime originalTimes(whichOrig,2)];
        
        % increase the block and start again
        whichOrig = whichOrig + 1;
        if whichOrig == size(originalTimes,1) + 1
            break
        end
        
        startTime = originalTimes(whichOrig,1);
        whichExc = 1;
    
    end
end

% Do a check for intersections
intersection = 0;
for i = 1:size(newTimes,1)
    for j = 1:size(excludedTimes,1)
        if doTimeRangesIntersectForgiving(newTimes(i,:),...
                excludedTimes(j,:),foregiveness) == 1
            intersection = 1;
        end
    end
    
end

if intersection == 1
    error('What\n');
end

% Remove any rows where the first column is bigger than the 2nd
removeRows = [];
for i = 1:size(newTimes,1)
    if newTimes(i,1) > newTimes(i,2)
        removeRows = [removeRows,i];
    end
end

newTimes(removeRows,:) = [];

end