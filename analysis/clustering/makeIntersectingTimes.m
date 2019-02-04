function newTimes = makeIntersectingTimes(time1,time2)

newTimes = [];

for i = 1:size(time1,1)
    
    for j = 1:size(time2,1)
        if doTimeRangesIntersectForgiving(time1(i,:),time2(j,:),0) == 1
            
            newTimes = [newTimes;max(time1(i,1),time2(j,1)),...
                min(time1(i,2),time2(j,2))];
            
        end
    
    end
end


end