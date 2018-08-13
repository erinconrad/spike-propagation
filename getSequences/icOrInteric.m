function ictal = icOrInteric(sequences,szTimes)

nseq = size(sequences,2);
ictal = zeros(nseq,1);

for i = 1:nseq
    curr_seq = sequences(:,i);
    B = sort(curr_seq);
    first_time = B(1);
    if first_time >= szTimes(1) && first_time <= szTimes(2)
        ictal(i) = 1;
    end
        
    
end

end