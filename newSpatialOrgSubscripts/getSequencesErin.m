function [sequences,discarded] = ...
    getSequencesErin(spikes, xyChan, t1, t2, minSeqLength, uniquenessCheck,...
    minUniqueSeqLength, maxPercTies,...
    minConcurrentFreqAbs, minConcurrentFreqRel, maxDist, minSpikesCloseEnough)

%{
 
    Function: detect sequences in spikes. Sequences are defined as spike trains
    that 1) contain minSeqLength or more spikes, 2) have latencies within
    t1 seconds of lead spike OR t2 seconds of prev spike. 
          
%}

% ------------------------------------------------------------------------------

discarded.origNum = 0; % number of original sequences pre-rejection
discarded.total =  0; % total number of rejected sequences
discarded.length = 0; % number of sequences rejected because too short
discarded.ties = 0; % number of sequences rejected because too many ties
discarded.spatialRestriction = 0;

allSeqsTies = [];

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
    if spikes(row,2) - head(2) <= t1 || spikes(row,2) - currseq(end,2) <=t2   

            
        currseq = vertcat(currseq,spikes(row,:));
        

    % Terminating spike reached- match matrix lengths, add sequence of overall if >= minSeqLength steps
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
            percTies = n_ties/length(ticksCol);
            allSeqsTies = [allSeqsTies,percTies];
            
            % Uniqueness check
            % This is my criteria for potentially throwing out a sequence if there are
            % too many simultaneous spikes. If this is set to 0, I never throw any
            % sequences out. If this is set to 1, I use the minimum number of
            % non-simultaneous spikes criterion. If this is set to 2, I use the maximum
            % percentage that is ties criterion.
            acceptablyUnique = 0;
            if uniquenessCheck == 0
                acceptablyUnique = 1;
            elseif uniquenessCheck == 1
                if length(unique(ticksCol)) >= minUniqueSeqLength
                    acceptablyUnique = 1;
                end
            elseif uniquenessCheck == 2
                if n_ties < length(ticksCol)*maxPercTies
                    acceptablyUnique = 1;
                end
            end
            
           % If acceptably unique, accept the sequence
           if acceptablyUnique == 1
                
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


% Reorder tied spikes based on spatial location
overall = reorder_tiesErin(overall, xyChan);



% Partition sequences- apply spatial restrictions, store
part  = spatialConstraint(overall, xyChan, minConcurrentFreqAbs,...
    minConcurrentFreqRel, maxDist,...
    minSpikesCloseEnough);

% store sequences according to which is the lead channel
sequences = part;
discarded.spatialRestriction = discarded.spatialRestriction + size(overall,2)/2 - size(sequences,2)/2;


 discarded.total = discarded.total + discarded.spatialRestriction;
 discarded.totalPercTies = mean(allSeqsTies);


end