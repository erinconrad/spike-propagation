function p = testRates(num1,time1,num2,time2,nboot)



% How many total spikes
total_n = num1 + num2;

% How much total time
time1 = round(time1);
time2 = round(time2);
total_time = time1 + time2;

rate_diff = zeros(nboot,1);

% Bootstrap loops
for ib = 1:nboot
    
    
    
    n1_fake = 0;
    n2_fake = 0;
    
    % Make total_n random spikes, random times
    X = randi(total_time,total_n,1);
    
    % Get number of spikes in each time block
    n1_fake = sum(X<=time1);
    n2_fake = sum(X>time1);
    if n1_fake + n2_fake ~= total_n
        error('What\n');
    end
    
    % Compare spike rates
    rate_diff(ib) = n1_fake/time1 - n2_fake/time2;
    
end

% Get actual spike rate diff
rate_diff_real = num1/time1 -  num2/time2;

% find the number of fake spike rates that are more extreme
rate_diff = sort(rate_diff);

p = (1 + sum(abs(rate_diff)>abs(rate_diff_real)))/(1+nboot);


end