function [p,z] = poisson_mean_diff(count1,count2,time1,time2)

%% Method I found online
% https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/PASS/Tests_for_the_Difference_Between_Two_Poisson_Rates.pdf

lambda1 = count1/time1;
lambda2 = count2/time2;

z = (lambda1-lambda2)/sqrt(lambda1/count1 + lambda2/count2);
p = normcdf(z);

%% Now do bootstrap method
% Randomize 1 versus 2 
% HOW????


end