function numClustRemoved(cluster)

whichPts = [];

for i = 1:length(cluster)
    if cluster(i).k ~= length(cluster(i).bad_cluster)
        whichPts = [whichPts,i];
    end
end

removed = [];
remain = [];

for whichPt = whichPts
    remain = [remain;cluster(whichPt).k-length(cluster(whichPt).bad_cluster)];
    removed = [removed;length(cluster(whichPt).bad_cluster)];
end

fprintf('The mean number of remaining clusters was %1.1f (range %d-%d)\n',...
    mean(remain),min(remain),max(remain));

fprintf('The mean number of removed clusters was %1.1f (range %d-%d)\n',...
    mean(removed),min(removed),max(removed));



end