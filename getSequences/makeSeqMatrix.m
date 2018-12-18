function seq_matrix = makeSeqMatrix(sequences,nch,doMorph)

if doMorph == 1
    nseq = size(sequences,2)/4;
    seq_matrix = nan(nch,nseq);
    
    for i = 1:nseq
        ch_col = sequences(:,(i-1)*4 + 1);
        time_col = sequences(:,(i-1)*4 + 2);
        height_col = sequences(:,(i-1)*4 + 3);
        width_col = sequences(:,(i-1)*4 + 4);
        
        % Loop through the channels in the sequences
        for j = 1:size(sequences,1)

            % reached the last channel in the sequence
            if ch_col(j) == 0
                break
            end

            % if there is a channel listed, then fill the matrix for that
            % sequence i and that channel ch_col(j) with the time it is
            % activated in that sequence
            seq_matrix(ch_col(j),i) = time_col(j);
        end
        
    end
else
    nseq = size(sequences,2)/2;
    seq_matrix = nan(nch,nseq);

    %% Loop through the sequences
    for i = 1:nseq
        ch_col = sequences(:,(i-1)*2 + 1);
        time_col = sequences(:,i*2);

        % Loop through the channels in the sequences
        for j = 1:size(sequences,1)

            % reached the last channel in the sequence
            if ch_col(j) == 0
                break
            end

            % if there is a channel listed, then fill the matrix for that
            % sequence i and that channel ch_col(j) with the time it is
            % activated in that sequence
            seq_matrix(ch_col(j),i) = time_col(j);


        end

    end

end

end