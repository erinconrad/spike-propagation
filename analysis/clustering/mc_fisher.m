np = 100;
ntrial = 500;

fisher_p = zeros(ntrial, 1);
for itrial = 1:ntrial
    pfake = max(0,rand(np,1)-0.05);
    ff = -2*sum(log(pfake));
    fisher_p(itrial) = 1-chi2cdf(ff, 2*np);
end