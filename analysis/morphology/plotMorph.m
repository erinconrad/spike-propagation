function plotMorph(pt,whichPts)

for whichPt = whichPts
    
    szTimes = pt(whichPt).newSzTimes;
    
    % Get all spike channels, times, height, width
    seqs = pt(whichPt).data.sequences;
    nseq = size(seqs,2)/4;
    
    all_spikes = [];
    
    for i = 1:nseq
        curr_seq = seqs(:,(i-1)*4+1:(i-1)*4+4);
        curr_seq(curr_seq(:,1) == 0,:) = [];
        ch_col = curr_seq(:,1);
        time_col = curr_seq(:,2);
        height_col = curr_seq(:,3);
        width_col = curr_seq(:,4);
        
        all_spikes = [all_spikes;...
            ch_col,time_col,height_col,width_col];
        
    end
    
    figure
    plot(all_spikes(:,2)/3600,smooth(all_spikes(:,3),1))
    hold on
    for j = 1:size(szTimes,1) 
        yl = ylim;
        plot([szTimes(j,1) szTimes(j,1)]/3600,yl,'k','LineWidth',2);
    end
    
    figure
    plot(all_spikes(:,2)/3600,smooth(all_spikes(:,4),1))
    hold on
    for j = 1:size(szTimes,1) 
        yl = ylim;
        plot([szTimes(j,1) szTimes(j,1)]/3600,yl,'k','LineWidth',2);
    end
    
    
end


end