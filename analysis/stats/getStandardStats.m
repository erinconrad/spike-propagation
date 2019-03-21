function test_stat = getStandardStats(x,y,which_stat)

if strcmp(which_stat,'sr') == 1
    % Do sign rank stat
    
    abs_diff = abs(x-y);
    sign_diff = sign(x-y);
    
    % remove pairs with sign 0
    zero_diff = find(sign_diff == 0);
    abs_diff(zero_diff) = [];
    sign_diff(zero_diff) = [];
    
    % Order according to difference
    [abs_diff,I] = sort(abs_diff);
    sign_diff = sign_diff(I);
    
    % Rank, starting with smallest as 1
    [r, ~] = tiedrank(abs_diff);
    
    W = sum(sign_diff.*r);
    N = length(sign_diff);
    z = W/sqrt(N*(N+1)*(2*N+1)/6);
    p = normcdf(z);
    test_stat = W;
    
elseif strcmp(which_stat,'rs') == 1
    % Do rank sum stat
    
    % Group all observations and assign ranks
    belong = [ones(length(x),1);2*ones(length(y),1)];
    [r,~] = tiedrank([x;y]);
    
    
    % add up the ranks in x and the ranks in y
    rsum_x = sum(r(belong==1));
    rsum_y = sum(r(belong==2));
    n_x = length(x);
    n_y = length(y);

    % pick smaller sample
    if n_x<n_y
        smsample = x; 
        small_sum = rsum_x;
    else
        smsample = y;
        small_sum = rsum_y;
    end
    
    % Get U
    U = small_sum - length(smsample)*(length(smsample)+1)/2;
    
    
    test_stat = U;
   

end

end