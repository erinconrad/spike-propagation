function [sequences,discarded] = getSequencesErin(spikes, xyChan)

%{
 
    Function: detect sequences in spikes. Sequences are defined as spike trains
    that 1) contain 3 or more spikes, 2) have latencies within 10 ticks of
    lead spike OR 3 ticks of prev spike. This program performs spatial
    partitioning of the sequences and applies a modest noise-reduction
    procedure. 
          
%}

% ------------------------------------------------------------------------------

%% Parameters
t1 = .05; % max time from first spike (50 ms)
t2 = .015; % max time from preceding spike (15 ms)

% minimum sequence length (how many spikes per sequence)
minSeqLength = 3; 

% maximum percentage of a sequence that is ties (when the spike hits
% multiple channels at the same time
maxPercTies = 0.7; 

discarded.origNum = 0; % number of original sequences pre-rejection
discarded.total =  0; % total number of rejected sequences
discarded.length = 0; % number of sequences rejected because too short
discarded.ties = 0; % number of sequences rejected because too many ties
discarded.remaining = 0; % number of sequences remaining at the end

%% Make sequences
nchans = length(xyChan);

% Initialize segment variables
sequences   = cell(1,nchans); 
   

head     = spikes(1,:); % the first spike channel and time
currseq  = head; % start the current sequence with the first spike

% this will be an nx(2*m) array, where n is the number of spikes in a
% sequence (padded with 0s so as to have uniform length per sequence, m is
% the number of sequences, and then within each m, the first column has the
% spike channel and the second column has the spike time
overall  = [];

% Start at second row of spike and loop through all the spikes
for row = 2:size(spikes,1)

    % If the time of the next spike is less than t1 from the first spike or
    % if the time of the next spike is less than t2 from the preceding
    % spike, add this next spike to the sequence
    if spikes(row,2) - head(2) <= t1 | spikes(row,2) - currseq(end,2) <=t2   

            
        currseq = vertcat(currseq,spikes(row,:));
        

    % Terminating spike reached- match matrix lengths, add sequence of overall if >= 3 steps
    else
        
        % add it to total num of original sequences
        discarded.origNum = discarded.origNum + 1;
        
        % If it's a long enough spike sequence
        if size(currseq,1) >= minSeqLength

           
            % a column with the spike times in the current sequence
            ticksCol = currseq(:,2);
            
            ticksCol(ticksCol==0) = []; % remove zeros (shouldnt be any)
            
            % n_ties is the number of non-unique spike times in the
            % sequence (so if a spike hit 2 channels at the same time, or
            % ties)
            n_ties   = length(ticksCol) - length(unique(ticksCol));
            
            % require that the percentage of ties be less than maxPercTies
            % of the total sequence
            if n_ties < length(ticksCol)*maxPercTies  % else, toss sequence

                % Add dummy zeros, concatenate
                
                % the difference between the length of the current sequence
                % and the length of the overall matrix of sequences
                diff = length(currseq) - size(overall,1);
                
                % if the current sequence is shorter, pad it with 0s
                if diff < 0
                    currseq = vertcat(currseq, zeros(abs(diff),2));
                    
                % if the current sequence is longer, pad the overall matrix
                % with zeros
                elseif diff > 0
                    overall = vertcat(overall, zeros(diff,size(overall,2)));
                end
                
                % horizontally concatenate the matrix of overall sequences
                % and the new sequence
                overall = [overall, currseq];

            else
                discarded.ties = discarded.ties + 1;
                discarded.total = discarded.total + 1;
            
            end
        else
            discarded.length = discarded.length + 1;
            discarded.total = discarded.total + 1;
        end

        % Update head and currseq
        head = spikes(row,:);
        currseq = head;

    end

end

% ERIN MADE IT TO HERE

% Reorder tied spikes based on spatial location
% CHECK ON THIS
overall = reorder_tiesErin(overall, xyChan);

% Extract sequences for each chan, apply spatial partition restrictions
discarded.remaining = 0;
 for lead = 1:nchans   
 % loop through each channel

    % Identify indices of sequences led by channel
    allseqs = []; 
    hits    = find(overall(1,:) == lead);
    hits(mod(hits,2)==0) = []; % remove even elements (ie, time columns)

    % Collect channel AND ticks columns into temp matrix
    % so temp is an nx(m*2) matrix where n is the number of spikes in a
    % sequence (padded with zeros) and m is the number of sequences led by
    % this channel and then for each m, column 1 is the channel and column
    % 2 is the time
    temp = [];
    for col = hits
        temp = [temp, overall(:,col), overall(:,col+1)];
    end

    % Partition sequences- apply spatial restrictions, store
    part  = spatialConstraint(temp, xyChan);
    sequences{lead} = part;
    
    discarded.remaining = discarded.remaining + size(part,2)/2;

 end
    
 discarded.percentdiscarded = discarded.total/discarded.origNum;


end