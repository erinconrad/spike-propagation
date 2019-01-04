function allSubsets = getContinuousSubset(n_total,n_subset)

%{
This function gets a subset equal to n_subset CONTINUOUS hours from a
larger sample of n_total hours.

One issue is that I am pre-definining the hour-long bins rather than
starting with an arbitrary time and then looking n_subset hours after that
time. I don't think this will make a huge difference.
%}

% E.g. if subset is 60 hours and total dataset is 70 hours, then there are
% 11 subsets possible (ending with hour 60, 61, 62, 63, ... 70).
num_of_subsets = n_total - n_subset + 1;

allSubsets = zeros(num_of_subsets,n_subset);

for i = 1:num_of_subsets
    
    % e.g., first is 1:60, second is 2:61, last is 11:70
    allSubsets(i,:) = i:i+n_subset-1;
end


end