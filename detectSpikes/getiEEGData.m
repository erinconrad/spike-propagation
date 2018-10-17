function data = getiEEGData(dataName,channels,indices,pwname)

% This is a tool to return information from a specified iEEG dataset and a
% certain time and channel arrangement

% If this is called with channels = 0, then it will just return some basic
% information about the data set

%% Unchanging parameters
loginname = 'erinconr';
%pwname = 'eri_ieeglogin.bin';

n = 0;


%% Open and get data
% I am putting this whole thing in a try-catch block inside a while loop
% because sometimes the ieeg.org server fails to give data, and this would
% usually crash the program



while n == 0

    try
        session = IEEGSession(dataName, loginname, pwname);
        fs = session.data.sampleRate;
        channelLabels = session.data.channelLabels;

        if sum(channels) ~= 0
            times = indices/fs;
            values = session.data.getvalues(indices, channels);
        end
        
        n = 1;
        
    catch
       fprintf('Failed to retrieve ieeg.org data, trying again...\n'); 
       n = 0; 
    end


end


%% Create struct
if sum(channels) ~= 0 
    data.values = values;
    data.times = times;
end
data.fs = fs;
data.name = dataName;
data.chLabels = channelLabels;

session.delete;
clearvars -except data

end