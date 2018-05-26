%{
This function loops through desired patients and loops through various
tmuls and absthresh's and runs portVisualizeSpikes to see how the spike
detector works for that time period

%}


function spikeVerificationNew(P)

tmuls_to_try = [11 13 15];
absthresh_to_try = [200 300 400];

for i = 1:1%length(P)
    
    if isempty(P(i).ieeg_name) == 1
        fprintf('Missing ieeg_name for patient %s, skipping\n',P(i).name);
        continue
    end
   
    [startTimes,chnames] = spikes_by_eye(P(i).name);
    
    if isempty(startTimes) || isempty(chnames)
        fprintf('Missing start times or channel names for patient %s\n',P(i).name);
        continue
    end
    
    chIds = zeros(size(chnames));
    for ich = 1:length(chnames)
        [Lia,chIds(ich)] = ismember(chnames{ich},P(i).electrodeData.unignoredChs);
        if Lia == 0
            fprintf('Warning, could not find channel %s in the list of unignored channels for patient %s\n',...
                chnames{ich},P(i).name);
            error('');
        end
    end
       
    for k = tmuls_to_try

        for m = absthresh_to_try

            portSpikesMultiple(P,i,chIds,k,m,startTimes);

        end



    end
        
    
    
    
    
    
end



end