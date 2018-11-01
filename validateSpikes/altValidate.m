function altValidate(pt,whichPt,whichDetector,tmul,absthresh,times,chnames)



if isempty(chnames) == 0
chs = zeros(size(chnames));
for ich = 1:length(chnames)
    [Lia,chs(ich)] = ismember(chnames{ich},pt(whichPt).electrodeData.unignoredChs);
    if Lia == 0
        fprintf('Warning, could not find channel %s in the list of unignored channels for patient %s\n',...
            chnames{ich},pt(whichPt).name);
        error('');
    end
end
else
    chs = 1:length(pt(whichPt).channels);
end

thresh.tmul = tmul;
thresh.absthresh = absthresh;

[gdf,extraOutput] = getSpikesSimple(pt,whichPt,times,whichDetector,thresh,0);
values = extraOutput.values;

seconds = 15;
window = seconds*pt(whichPt).fs;
n_chunks = ceil((times(2)-times(1))/seconds);
%% Make plots
for i = 1:n_chunks
    indices = round([(i-1)*window+1:min(i*window,size(values,1))]);
    min_time = [times(1)+(i-1)*seconds];
    max_time = times(1) + min(i*seconds,size(values,1)/pt(whichPt).fs);
    pl_values = values(indices,:);
    pl_spikes = gdf(gdf(:,2)>=min_time & gdf(:,2)<= max_time,:);
    
    figure
    set(gcf,'Position',[200 200 1400 800])
    toAdd = 0;
   
    for ch = chs
        plot(linspace(min_time,max_time,length(indices)),...
            pl_values(:,ch)+repmat(toAdd,length(pl_values(:,ch)),1));
        hold on
        
        ch_pl_spikes = pl_spikes(pl_spikes(:,1) == ch,2);
        value_at_spike = values(round((ch_pl_spikes - times(1))*pt(whichPt).fs),ch);
        scatter(ch_pl_spikes,value_at_spike + repmat(toAdd,length(value_at_spike),1),60);
        
        if ch ~= chs(end)
            toAdd = toAdd + max(pl_values(:,ch)) - min(pl_values(:,ch+1))+ 30;
        end
        
    end
    
    fprintf('Press any key\n');
    pause;
    close(gcf)
    
    
end



end