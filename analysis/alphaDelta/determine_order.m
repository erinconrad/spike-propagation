function [p,t,b] = determine_order(X,Y,plotInfo)

%{
If no significant lags in the first 10, I don't do an AR model. If there
are any significant lags, I do an AR(2) model. I did this because when I
plotted the partial correlations, AR(2) was enough for most of the
patients.
%}

%% Do linear regression (no AR)
[b1,~,resid,~,stats] = regress(Y, X);

% Get initial guess for autocorrelation term
r = corr(resid(1:end-1),resid(2:end));  

%% Calculate the lags with significant partial correlation (order of the AR model)
[pacf,lags,bounds] = parcorr(resid,'NumLags',10);
sig_lags = lags(find(abs(pacf) > bounds(1)));
sig_lags(sig_lags == 0) = [];
n_lags = length(sig_lags);

if n_lags > 0
    n_lags = 2;
end

%% Plot the autocorrelation and partial correlation of the residuals
if plotInfo == 1
    figure
    parcorr(resid,'NumLags',10)
    
    fprintf('AR lags are %d.\n',n_lags);
    pause
    close(gcf)
end

%% Do appropriate model based on number of lags
if n_lags == 0
    % Not an AR model
    
    % Do linear regression
    mdl = fitlm(X,Y);
    t = mdl.Coefficients.tStat(2);
    p = mdl.Coefficients.pValue(2);
    b = mdl.Coefficients.Estimate(2);
    constantB = mdl.Coefficients.Estimate(1);
    new_LL = mdl.LogLikelihood;
    num_pred = mdl.NumPredictors;
    new_AIC = 2*num_pred - 2*new_LL;
else
    % AR model of order n_lags
    
    % Make ARIMA model
    n_lags = 2;
    mdl = regARIMA('ARLags',[1,2]);
    %mdl = regARIMA('ARLags',sig_lags');
    mdl = estimate(mdl,Y,'beta0',b1(1),'Intercept0',b1(2),...
    'X',X(:,1),'Display','Off');
    %mdl = regARIMA('ARLags',n_lags);
    %{
    mdl = estimate(mdl,Y,'beta0',b1(1),'Intercept0',b1(2),...
    'X',X(:,1),'AR0',r,'Display','Off');
    %}
    results = summarize(mdl);
    p = results.Table.PValue(2+n_lags);
    t = results.Table.TStatistic(2+n_lags);
    b = results.Table.Value(2+n_lags);
    constantB = results.Table.Value(1);
    new_LL = results.LogLikelihood;
    new_AIC = results.AIC;
    
    
end

%% Plot fits
if plotInfo == 1
    figure
    subplot(1,2,1)
    plot(Y,'r')
    hold on
    plot(X*b1,'k--');
    legend('Real Y','Model')
    title('Linear regression');
    
    subplot(1,2,2)
    plot(Y,'r')
    hold on
    plot(X*[b;constantB],'k--')
    legend('Real Y','Model')
    title('Autoregressive model (if appropriate')
    
    % Get AIC
    old_mdl = fitlm(X,Y);
    old_LL = old_mdl.LogLikelihood;
    num_pred = old_mdl.NumPredictors;
    old_AIC = 2*num_pred - 2*old_LL;
    
    LL_rat = new_LL-old_LL; %log(A/B) = logA ? logB so this is log(LR_new/LR_old)
    if LL_rat > 0
        fprintf('LL favors new model:\n old LL is %1.2f and new LL is %1.2f\n',...
            old_LL,new_LL);
    elseif LL_rat < 0
        fprintf('LL favors old model:\n old LL is %1.2f and new LL is %1.2f\n',...
            old_LL,new_LL);
    else
        fprintf('LL favors nothing:\n old LL is %1.2f and new LL is %1.2f\n',...
            old_LL,new_LL);
    end
    
    if new_AIC < old_AIC
        fprintf('AIC favors new model:\n old AIC is %1.2f and new AIC is %1.2f\n',...
            old_AIC,new_AIC);
    elseif old_AIC < new_AIC
        fprintf('AIC favors old model:\n old AIC is %1.2f and new AIC is %1.2f\n',...
            old_AIC,new_AIC);
    else
        fprintf('AIC favors nothing:\n old AIC is %1.2f and new AIC is %1.2f\n',...
            old_AIC,new_AIC);
    end
    
    fprintf('Old b is %1.1f and new b is %1.1f\n',b1(1),b);
    fprintf('Old p is %1.2e and new p is %1.2e\n',stats(3),p);
    
    pause
    close(gcf)
end

% Old
%{
MdlY = arima('AR',results.Table.Value(2),'Constant',results.Table.Value(1),...
    'Variance',results.Table.Value(4),'Beta',[results.Table.Value(3)]);

Y_fake = simulate(MdlY,size(X,1),'X',X(:,1));


%% Define anonymous function representing new fit including autocorrelation
f = @(c,x) [Y(1); c(1)*Y(1:end-1) + (x(2:end,:)- c(1)*x(1:end-1,:))*c(2:end)];

mdl = fitnlm(X,Y,f,[r;b]);
p_old = mdl.Coefficients.pValue(2);
%r2 = mdl.Rsquared.Adjusted;
t_old = mdl.Coefficients.tStat(2);
b_old = mdl.Coefficients.Estimate(2);



if plotInfo == 1
    figure
    subplot(1,3,1)
    plot(Y,'r')
    hold on
    plot(X*[b_old;mdl.Coefficients.Estimate(3)],'k--')
    title('Old AR(1) model')

    subplot(1,3,2)
    plot(Y,'r')
    hold on
    plot(X*[results.Table.Value(3);results.Table.Value(1)],'k--');
    title('New AR(1) model')
    
    subplot(1,3,3)
    plot(Y,'r')
    hold on
    plot(EstMdl,'k--');
    title('New AR(1) model via simulate')
    
    
    % Show the old and new coefficients
    [b_old results.Table.Value(3);
     mdl.Coefficients.Estimate(3) results.Table.Value(1)]
    



end

%}

%{
https://newonlinecourses.science.psu.edu/stat462/node/189/

- Partial autocorrelation function assess the appropriate lags for the
errors in a regression model with autoregressive errors. Fit a regression
model to the time series data and plot the partial autocorrelation for the
residuals versus the lag.



%}

%{
% Do linear regression
predictor = [ones(size(X)) X];
%[b,~,resid,~,stats] = regress(Y, X);
beta = predictor\Y;
u = Y - predictor*beta;

figure
plot(u)
pause
close(gcf)

figure
subplot(2,1,1)
autocorr(u,'NumLags',200)
subplot(2,1,2)
parcorr(u)
pause
close(gcf)

dlY = diff(log(Y));
dlX = diff(log(X));
dlPredictor = [ones(size(dlX)) dlX];
beta = dlPredictor\dlY;
u = dlY - dlPredictor*beta;
figure
plot(u);
pause
close(gcf)

figure
subplot(2,1,1)
autocorr(u)
subplot(2,1,2)
parcorr(u)
pause
close(gcf)
%}



end