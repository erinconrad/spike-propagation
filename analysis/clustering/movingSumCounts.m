function [counts,window_times] = movingSumCounts(times,all_times,window)


% get min and max times
totalTime = [floor(min(all_times)),ceil(max(all_times))];

% initialize counts and times
if iscell(times) == 1
    counts = zeros(length(times),totalTime(2) - totalTime(1) + 1);
else
    counts = zeros(1,totalTime(2) - totalTime(1) + 1);
end
window_times = zeros(1,totalTime(2) - totalTime(1) + 1);


% Loop through all times
n = 0;
for t = totalTime(1):totalTime(2)
    n = n+1;
    % Time range goes from t to t + window
    time_range = [t,min(t + window,totalTime(2))];
    
    
    %window_times(n) = (time_range(1)+time_range(2))/2;
    
    % The value in just reporting the last time is that if I want to see if
    % something changes before a seizure, this will let me know the last
    % time point I am averaging to see if I am including any post-seizure
    % business
    window_times(n) = time_range(2);

    % Loop through clusters
    if iscell(times) == 1
    
        for j = 1:length(times)
           counts(j,n) = sum(times{j}>=time_range(1) & ...
               times{j}<= time_range(2));
        end
    else
       counts(n) = sum(times>=time_range(1) & times<=time_range(2)); 
    end

end

end