function intersect = doTimeRangesIntersect(a,b)

% Tests if 2 time ranges a and b intersect

% By default, they intersect
intersect = 1;

%% Ways to not intersect

% a(2) is less than b(1)
if a(2) < b(1)
    intersect = 0;
end

% if a(1) is bigger than b(2)
if a(1) > b(2)
    intersect = 0;
end


end