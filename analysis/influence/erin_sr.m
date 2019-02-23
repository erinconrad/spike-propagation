function erin_sr(x,y)

if size(x,2) > size(x,1), x=x'; end
if size(y,2) > size(y,1), y=y'; end

abs_diff = abs(y-x);
sign_diff = sign(y-x);

abs_diff(sign_diff == 0) = [];
sign_diff(sign_diff == 0) = [];

[~,I] = sort(abs_diff);

N = length(abs_diff);
W = sum(sign_diff(I).*[1:N]');
ranked = tiedrank(abs_diff);
sum_pos = sum(ranked(sign_diff > 0));
sum_neg = sum(ranked(sign_diff < 0));
T = min(sum_pos,sum_neg);


if N >= 20
    sigma = sqrt(N*(N+1)*(2*N+1)/6);
    z = W/sigma;
    p = 2*normcdf(-abs(z));
end

end