function data =  alphaDelta(pt,whichPts)

%% Load file paths, etc.
[~,~,~,~,pwfile] = fileLocations;

for whichPt = whichPts
    
    fs = pt(whichPt).fs;
    dataName = pt(whichPt).ieeg_name;
    channels = pt(whichPt).channels;
    nch = length(channels);
    
    %  Loop over run times
    for tt = 1:size(pt(whichPt).runTimes,1)
        
        % Get the desired indices
        desiredTimes = pt(whichPt).runTimes(tt,:);
        indices = round(desiredTimes(1)*fs):round(desiredTimes(2)/20*fs);
        
        % Get the data
        tic
        data = getiEEGData(dataName,channels,indices,pwfile);
        toc
        fprintf('Retrieved data\n');
        
        % remove nans
        data.values(isnan(data.values)) = 0;
        
        % Loop over channels
        %{
        for dd = 1:nch
            X = data.values(:,dd);
            
            % Calculate fft
            Y = fft(X);
            
            Y;
            
        end
        %}
        
       
    
    end
end


end