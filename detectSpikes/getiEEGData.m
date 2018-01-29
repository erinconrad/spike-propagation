function data = getiEEGData(dataName,channels,indices,pwname)

% This is a tool to return information from a specified iEEG dataset and a
% certain time and channel arrangement

% If this is called with channels = 0, then it will just return some basic
% information about the data set

%% Unchanging parameters
loginname = 'erinconr';
%pwname = 'eri_ieeglogin.bin';

%% Open and get data
session = IEEGSession(dataName, loginname, pwname);
fs = session.data.sampleRate;
channelLabels = session.data.channelLabels;

if sum(channels) ~= 0
    times = indices/fs;
    values = session.data.getvalues(indices, channels);
end


%% Create struct
if sum(channels) ~= 0 
    data.values = values;
    data.times = times;
end
data.fs = fs;
data.name = dataName;
data.chLabels = channelLabels;

end