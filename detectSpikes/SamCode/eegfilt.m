function out = eegfilt(x,fc,type, fs)
%filter eeg data using Butterworth filter
%out = eegfilt(data,cutfreq,typ);
%out = eegfilt(data,70,'hp'); high pass with 70Hz cutoff

%EEG_BUTTER - Butterworth filter implementation
%  xf = eeg_butter(x,sampl_freq,cutoff_freq,filter_type,num_poles)

np = 6;
if sum(fc >= fs/2), error('Cutoff frequency must be < one half the sampling rate'); end
fn = fs/2;
type = type(1:2);
if strcmp(type,'bp'), type = 'lp'; end

switch type,
case 'lp',
   [B,A] = butter(np,fc/fn);
case 'hp',
   [B,A] = butter(np,fc/fn,'high');
case 'st'
   [B,A] = butter(np,fc/fn,'stop');
end

out = filtfilt(B,A,x);
