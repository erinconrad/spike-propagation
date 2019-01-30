function SD = standardDistance(locs)

%{

Calculate the standard distance between a bunch of points (a measure of
spatial dispersion that should not depend on the number of points

%}

n = size(locs,1);

SD = sqrt((...
    sum((locs(:,1)-mean(locs(:,1))).^2) + ...
    sum((locs(:,2)-mean(locs(:,2))).^2) + ...
    sum((locs(:,3)-mean(locs(:,3))).^2))/...
    n);


end