function [r,p,h,ci] = pointBiserial(X,Y,alpha)

%{
Point biserial correlation coefficient

X is continuous (recruitment latency)
Y is binary (SOZ or not)

Divide data set into 2 groups, group 1 which received the value "1" on Y,
and group 2, which received the value "0" on Y.


%}

n1 = sum(Y==1); % number in group 1
n0 = sum(Y==0); % number in group 2
n = n0 + n1; % total number

M1 = mean(X(Y==1)); % mean of X in group 1
M0 = mean(X(Y==0)); % mean of X in group 2

s = sqrt(1/n*sum((X-mean(X)).^2)); % standard deviation of whole population of X

r = (M1-M0)/s*sqrt(n1*n0/n^2); % Point biserial correlation coefficient

[p,h,ci] = ranksum(X(Y==1),X(Y==0),alpha); % non-parametric Wilcoxon rank sum test

end