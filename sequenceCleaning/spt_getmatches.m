function [ matches,dist_matches ] = spt_getmatches( point,ref_seq,tst_seq,ss_thresh,tt_thresh )

    % Get temporal matches in reference sequence
    point_time   = tst_seq(point,4);
    offsets      = abs(ref_seq(:,4) - point_time);
    temp_matches = find(offsets <= tt_thresh);
    
    % Get spatial matches in reference sequence
    distances    = spt_dist(tst_seq(point,1:3)',ref_seq(:,1:3)');
    spat_matches = find(distances(1,:)<=ss_thresh);
    
    % Get overlap, return distances for matching chans
    matches      = intersect(temp_matches,spat_matches);
    dist_matches = distances(1,matches);
    
end

