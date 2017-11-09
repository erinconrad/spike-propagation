function Patient = getRecruitmentLatency(Patient)

nSegments = size(Patient.sequences,2); % Number of segments
nChannels = size(Patient.xyChan,1); % Number of channels

recruitmentLatencySingle = cell(1,nSegments);
spikeCount = cell(1,nSegments);


for iSeg = 1:nSegments
    recruitmentLatencySingle{iSeg} = cell(1,nChannels);
    spikeCount{iSeg} = zeros(1,nChannels);
end

for iSeg = 1:nSegments
   for iChannel = 1:nChannels
       
       
       % Initialize recruitmentLatencySingle, which for each segment and
       % each channel, will be an nx2 array where n is the number of total
       % spike sequences
       recruitmentLatencySingle{iSeg}{iChannel} = nan(length(Patient.sequences{iSeg}),2);  
   end
end



for iSeg = 1:nSegments
    for iChannel = 1:nChannels
        if isempty(Patient.sequences{iSeg}{1,iChannel}) == 1 % If there are no sequences for that segment and channel, continue
            continue
        else
            nSeqs = size(Patient.sequences{iSeg}{1,iChannel},2)/2; % The number of sequences is the number of columns divided by 2
            
            % Loop through all the sequences
            for jSeq = 1:nSeqs
                column = (jSeq-1)*2+1; % This is the column showing the channel
                
                headChannel = Patient.sequences{iSeg}{1,iChannel}(1,column); % The first row of the spike sequence contains the first spike
                headTime = Patient.sequences{iSeg}{1,iChannel}(1,column+1); % The next column is the spike time
                
                for kSpikeInSeq = 1:size(Patient.sequences{iSeg}{1,iChannel},1) % The number of rows in the sequence (number of spikes, padded with zeros)
                    
                    % If the spike data is zero, continue
                    if Patient.sequences{iSeg}{1,iChannel}(kSpikeInSeq,column) == 0
                        continue
                    else
                        
                        % Get the latency at which that channel was
                        % activated relative to the head channel
                        tempLatency = Patient.sequences{iSeg}{1,iChannel}(kSpikeInSeq,column+1) - headTime;
                        
                        % Get the channel being activated
                        tempChan = Patient.sequences{iSeg}{1,iChannel}(kSpikeInSeq,column);
                        
                        % Increase the spike count for that channel
                        spikeCount{iSeg}(tempChan) = spikeCount{iSeg}(tempChan) + 1;
                        tempSCount = spikeCount{iSeg}(tempChan);
                        
                        % Fill up the first column with the latency
                        recruitmentLatencySingle{iSeg}{tempChan}(tempSCount,1) = tempLatency; 
                        
                        % Fill up the second column with the head channel
                        recruitmentLatencySingle{iSeg}{tempChan}(tempSCount,2) = headChannel;
                    end
                end
            end
        end
    end
end

%% Trim NaNs
for iSeg = 1:nSegments
    for iChannel = 1:nChannels
        [rownan,~] = find(isnan(recruitmentLatencySingle{iSeg}{iChannel}));
        
        
        % Do nothing if there are no nans
        if isempty(rownan) == 1
            
            
        % if the first row has nans then they're all nans so there are no
        % spikes in that channel, so replace it with an empty matrix
        elseif rownan == 1
            recruitmentLatencySingle{iSeg}{iChannel} = [];
            
        % Otherwise discard all the non nans
        else 
            recruitmentLatencySingle{iSeg}{iChannel} = recruitmentLatencySingle{iSeg}{iChannel}(1:rownan-1,:);
            
        end
            
    end
    
end


Patient.recruitmentLatencySingle = recruitmentLatencySingle;
Patient.spikeCount = spikeCount;




end