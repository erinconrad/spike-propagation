function [clueGraph] = ...
    spt_clueaware(ss_thresh,tt_thresh,xlocs,ylocs,zlocs,chanlocs,tstamps)

% This function uses the 'Clue-Aware' trajectory clustering algorithm
% presented in Hung et al. 2011 (doi: 10.1007/s00778-011-0262-6). The
% result is a pairwise similarity measure for spike sequences based on
% spatial and temporal thresholding values (set in the main script)


clueGraph   = zeros(size(chanlocs,2),size(chanlocs,2));
for ref_idx = 1:size(xlocs,2)
    fprintf('%1.2f percent complete\n',ref_idx/size(xlocs,2)*100);
       
    
    % Initialize reference sequence and similarity score to test sequence
    ref_seq = spt_initseq(xlocs,ylocs,zlocs,tstamps,ref_idx);
    
    for tst_idx = 1:size(xlocs,2)

        totalScore        = 0; 
        tst_seq           = spt_initseq(xlocs,ylocs,zlocs,tstamps,tst_idx);       
        
        for point = 1:size(tst_seq,1)
            
            % Find spatial and temporal matches in reference seq
            [matches,dist_matches] = spt_getmatches(point,ref_seq,tst_seq,ss_thresh,tt_thresh);
            
            if ~isempty(matches)
            
                % Compute Similarity Score for all matches
                simScores  = 1 - dist_matches./ss_thresh;
                totalScore = totalScore + max(simScores);
                
            end
            
        end
        
        % Store sequence similarity in clueGraph
        % NOTE: row = reference sequence, col = test sequence
        clueGraph(ref_idx,tst_idx) = totalScore/length(tst_seq);
        
    end    
end         
end                
    
    
    
    
    
    
    
    