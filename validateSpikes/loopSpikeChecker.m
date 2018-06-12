function pt=loopSpikeChecker(pt)

tmuls_to_try = [15];
absthresh_to_try = [300];

for i = 1%:length(pt)
    
    if isempty(pt(i).ieeg_name) == 1
        fprintf('Missing ieeg_name for patient %s, skipping\n',pt(i).name);
        continue
    end
   
    %[startTimes,chnames] = spikes_by_eye(pt(i).name);
    startTimes = [10801.11; 10861.95; 10880.03; 10881.60;42138.89;...
         42074.98; 42076.13; 42120.25;42208.89; 42347.86;...
          42626.85; 42647.54;42847.97; 43025.92;86426.04;...
         86841.83;87300.85;87318.10;87460.41;87461.30];
     chnames = {'LG2','LG3','LG4','LG10','LG17','LG25',...
        'LG34','LG52','LG54','LG59','LIH5','RFR3','LST1',...
        'LST2','LST3'};
    
    if isempty(startTimes) || isempty(chnames)
        fprintf('Missing start times or channel names for patient %s\n',P(i).name);
        continue
    end
    
    chIds = zeros(size(chnames));
    for ich = 1:length(chnames)
        [Lia,chIds(ich)] = ismember(chnames{ich},pt(i).electrodeData.unignoredChs);
        if Lia == 0
            fprintf('Warning, could not find channel %s in the list of unignored channels for patient %s\n',...
                chnames{ich},pt(i).name);
            error('');
        end
    end
       
    for k = tmuls_to_try

        for m = absthresh_to_try

            [pt(i).sensitivity,pt(i).accuracy] = spikeChecker(pt,i,chIds,...
    startTimes,ones(length(startTimes),1),k,m);

        end



    end
        
    
    
    
    
    
end



end


