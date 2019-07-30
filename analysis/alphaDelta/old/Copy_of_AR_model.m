function [p,t,b_out] = AR_model(X,Y,plotInfo)

% Do linear regression
[b,~,resid,~,stats] = regress(Y, X);

p_AR = dwtest(resid,X);

% Get initial guess for autocorrelation term
r = corr(resid(1:end-1),resid(2:end));  

% Define anonymous function representing new fit including autocorrelation
f = @(c,x) [Y(1); c(1)*Y(1:end-1) + (x(2:end,:)- c(1)*x(1:end-1,:))*c(2:end)];
[c,~,~,CovB,~,~] = nlinfit(X,Y,f,[r;b]);
mdl = fitnlm(X,Y,f,[r;b]);
p = mdl.Coefficients.pValue(2);
%r2 = mdl.Rsquared.Adjusted;
t = mdl.Coefficients.tStat(2);
b_out = c(2);

%{
- autocorr
- AR model vs MA model

- reviewer 
- Tom Nichols at Oxford
- Tim Johnson at Michigan in biostats
- Hernando Ombao at Kaust

%}

if plotInfo == 1     
    
    fprintf('P-value for DW test: %1.1e.\n',p_AR);

    %Info about new residuals
    figure
    subplot(1,3,1)
    u = Y - f(c,X);
    plot(u);
    title('Residuals');

    subplot(1,3,2)
    plot(Y,'b');
    hold on
    plot(X*b,'r--');
    legend('Real thing','original model');

    subplot(1,3,3)
    plot(Y,'b');
    hold on
    plot(f(c,X),'r--');
    legend('Real','AR model');
    pause
    close(gcf)
end

end