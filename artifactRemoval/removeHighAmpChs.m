function remove = removeHighAmpChs(data,fs)
%% parameters
fr = 2;


%% Low pass filter the data
ldata = eegfilt(data,fr,'lp',fs);


end