function time_chunks = getUnremovedTimes(pt,whichPt)

% NEED TO INSERT A REMOVED ARRAY FOR EACH ALLTIME
    

time_chunks = [];

% Loop through all run times
for i = 1:size(pt(whichPt).allTimes,1)
    if isempty(removed) == 1
        time_chunks = [time_chunks;pt(whichPt).allTimes(i,:)];
    else        

        % Add time from beginning of run to first removed
        time_chunks = [time_chunks;...
            pt(whichPt).allTimes(i,1) removed(1,1)];

        if size(removed,1) > 1

            % then add times from end of last removed to beginning of
            % next removed
            for j = 2:size(removed,1)
                time_chunks = [time_chunks;...
                    removed(j-1,2) removed(j,1)];

            end

            % then add the remainder going to the end of the run
            time_chunks = [time_chunks;...
                removed(size(removed,1),2) pt(whichPt).allTimes(i,2)];
        else
            time_chunks = [time_chunks;...
                removed(1,2) pt(whichPt).allTimes(i,2)];
        end

    end

end



end