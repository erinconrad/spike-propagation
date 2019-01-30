function testAR

%% make random data
X = rand(100,1);
Y = rand(100,1);

%% Add a cyclical trend
slope = sin(1:100)';

X = X + slope;
Y = Y + slope;

%% Plot them
figure
subplot(2,1,1)
plot(X)

subplot(2,1,2)
plot(Y)

%% Do initial regression
X = [ones(size(X)) X];
[b,~,resid,~,stats] = regress(Y, X);

fprintf('The p value for linear regression is %1.1e\n',stats(3));

%% Get initial guess for autocorrelation term
r = corr(resid(1:end-1),resid(2:end));  
    
 %% Define anonymous function representing new fit including autocorrelation
f = @(c,x) [Y(1); c(1)*Y(1:end-1) + (x(2:end,:)- c(1)*x(1:end-1,:))*c(2:end)];
mdl = fitnlm(X,Y,f,[r;b]);

fprintf('The new p value is %1.1e\n',mdl.Coefficients.pValue(3))

    

end