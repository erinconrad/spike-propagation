function ad_rat = innerAlphaDelta(dataName,channels,indices,pwfile,indicesToClip,fs,whichVer)

%{
Calculates the alpha delta ratio
%}

% Get the data
tic
data = getiEEGData(dataName,channels,indices,pwfile);
%pause(10)
toc
fprintf('Retrieved data, doing analysis.\n');


% remove nans
data.values(isnan(data.values)) = 0;

% Remove seizure times
data.values(max(indicesToClip-indices(1),1),:) = [];
nch = length(channels);

ad_rat = zeros(nch,1);
% Loop over channels

for dd = 1:nch
    X = data.values(:,dd);

    % subtract mean
    X = X - mean(X);


    %% fft approach

    % Calculate fft
    Y = fft(X);

    % Get power
    P = abs(Y).^2;
    freqs = linspace(0,fs,length(P)+1);
    freqs = freqs(1:end-1);

    % Take first half
    P = P(1:ceil(length(P)/2));
    freqs = freqs(1:ceil(length(freqs)/2));

    %plot(freqs,P);

    alpha = sum(P(freqs>=8 & freqs<=13));
    delta = sum(P(freqs>=1 & freqs<=4));

    % Get alpha/delta ratio
    if whichVer == 1
        ad_rat(dd) = alpha/delta;
    elseif whichVer == 2
        ad_rat(dd) = bandpower(X,fs,[8 13])/bandpower(X,fs,[1 4]);
    elseif whichVer == 3
        ad_rat(dd) = alpha;
    elseif whichVer == 4
        ad_rat(dd) = delta;
    end

    %fprintf('Got FFT approach\n');

    % Get all frequency bands
    %{
    for bb = 1:nbands
        frange = (bb-1)*10:min(bb*10,fs/2);
        all_p(dd,bb,tt) = sum(P(freqs>=frange(1) & freqs<=frange(2)));

    end
    %}

    %% Bandpass approach (gives about the same result as above, but slower)



   % ad_rat_band(dd,tt) = bandpower(X,fs,[8 13])/bandpower(X,fs,[1 4]);




end

clearvars -except ad_rat

end