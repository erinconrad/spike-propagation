function [p,z] = poisson_mean_diff(count1,count2,time1,time2)

%% Method I found online
% https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/PASS/Tests_for_the_Difference_Between_Two_Poisson_Rates.pdf

lambda1 = count1/time1;
lambda2 = count2/time2;

z = (lambda1-lambda2)/sqrt(lambda1/count1 + lambda2/count2);
p = normcdf(z);

true_diff = lambda1-lambda2;

%% Now do bootstrap method
% Randomize 1 versus 2 
nb = 1e3;
diff_boot = zeros(nb,1);
for ib = 1:nb
    c1 = randi(count1+count2);
    c2 = count1+count2 - c1;
    diff_boot(ib) = c1/time1-c2/time2;
    
end

sort_diff = sort(diff_boot);
%{
figure
plot(sort_diff)
%}

end