%{

This function takes a bunch of spike sequences and calculates recruitment
latency, which is how delayed a channel tends to be activated in the spike
sequence

%}

function [recruitmentLatencySingle,spikeCount] = getRecruitmentLatency(sequences,xyChan)

nChannels = size(xyChan,1); % Number of channels

% initialize output variable
recruitmentLatencySingle = cell(1,nChannels);
spikeCount = zeros(1,nChannels);


for iChannel = 1:nChannels


   % Initialize recruitmentLatencySingle, which for 
   % each channel, will be an nx2 array where n is the number of total
   % spike sequences. It will say when in that particular spike sequence
   % the channel was activated relative to the lead channel. This will
   % mostly be nans because most channels are not activated in a given
   % sequence.
   recruitmentLatencySingle{iChannel} = nan(length(sequences),2);  
end

nSeqs = size(sequences,2)/2; % The number of sequences is the number of columns divided by 2

% Loop through all the sequences
for jSeq = 1:nSeqs
    column = (jSeq-1)*2+1; % This is the column showing the channel

    headChannel = sequences(1,column); % The first row of the spike sequence contains the first spike
    headTime = sequences(1,column+1); % The next column is the spike time

    for kSpikeInSeq = 1:size(sequences,1) % The number of rows in the sequence (number of spikes, padded with zeros)
    % Loop through each spike in the sequence

        % If the spike data is zero, continue
        if sequences(kSpikeInSeq,column) == 0
            continue
        else

            % Get the latency at which that channel was
            % activated relative to the head channel
            tempLatency = sequences(kSpikeInSeq,column+1) - headTime;

            % Get the channel being activated
            tempChan = sequences(kSpikeInSeq,column);

            % Increase the spike count for that channel that is activated
            spikeCount(tempChan) = spikeCount(tempChan) + 1;
            tempSCount = spikeCount(tempChan);

            % Fill up the first column with the latency
            % This is basically just saying that for the channel of
            % interest (NOT THE LEAD CHANNEL), I want to look at
            % the next free row in the recruitmentLatencySingle
            % array for that channel and add the latency and the
            % head channel number for this particular appearance in
            % a sequence
            recruitmentLatencySingle{tempChan}(tempSCount,1) = tempLatency; 

            % Fill up the second column with the head channel
            recruitmentLatencySingle{tempChan}(tempSCount,2) = headChannel;
        end
    end
end




%% Trim NaNs

for iChannel = 1:nChannels
    [rownan,~] = find(isnan(recruitmentLatencySingle{iChannel}));


    % Do nothing if there are no nans
    if isempty(rownan) == 1


    % if the first row has nans then they're all nans so there are no
    % spikes in that channel, so replace it with an empty matrix
    elseif rownan == 1
        recruitmentLatencySingle{iChannel} = [];

    % Otherwise discard all the nans
    else 
        recruitmentLatencySingle{iChannel} = recruitmentLatencySingle{iChannel}(1:rownan-1,:);

    end

end
    



recruitmentLatencySingle = recruitmentLatencySingle;
spikeCount = spikeCount;




end