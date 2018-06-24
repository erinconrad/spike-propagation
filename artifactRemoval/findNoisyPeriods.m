%{
 This function calculates some measure of noise in the data and gives a
 binary 1 (if it's "noisy") or 0 (if it's not)


%}

function noise_bin = findNoisyPeriods(data,whichMethod)

noise_bin = -1;

%% diff method
if whichMethod == 1
    % Define the noise as the sqrt of the sum of the squared difference between
    % adjacent time points (so basically the amplitude of the high
    % frequency content of the data)
    noise = sqrt(sum(diff(data).^2));
    if noise > 1e4
        noise_bin = 1;
    end

%% RMS method
elseif whichMethod ==2
    % Define the noise as the RMS of the data
    noise = rms(data);
    if noise > 400
        noise_bin = 1;
    end
    
end


end