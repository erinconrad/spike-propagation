function ad_rat = innerAlphaDelta(dataName,channels,indices,pwfile,indicesToClip)

     % Get the data
    tic
    data = getiEEGData(dataName,channels,indices,pwfile);
    toc
    fprintf('Retrieved data, doing analysis.\n');


    % remove nans
    data.values(isnan(data.values)) = 0;

    % Remove seizure times
    data.values(indicesToClip-indices(1),:) = [];
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

        % Get alpha/delta ratio
        alpha = sum(P(freqs>=8 & freqs<=13));
        delta = sum(P(freqs>=1 & freqs<=4));
        ad_rat(dd) = alpha/delta;

        %fprintf('Got FFT approach\n');

        % Get all frequency bands
        %{
        for bb = 1:nbands
            frange = (bb-1)*10:min(bb*10,fs/2);
            all_p(dd,bb,tt) = sum(P(freqs>=frange(1) & freqs<=frange(2)));

        end
        %}

        %% Bandpass approach



       % ad_rat_band(dd,tt) = bandpower(X,fs,[8 13])/bandpower(X,fs,[1 4]);




    end

end