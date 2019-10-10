function power = alphaDelta(whichPts,pt)

%{
this will load a pt struct with
spike sequence info (although it will really only use the run times from
this struct) and it will grab ieeg data and calculate the alpha delta power
ratio in 2000 second bins). It calls innerAlphaDelta to do the actual work
(which is in the same folder). This will output a power struct with the
alpha delta power ratio.
%}

if exist('pt','var') == 0
    example = 0;
%% Load file paths, etc.
[~,~,scriptFolder,resultsFolder,pwfile] = fileLocations;
p1 = genpath(scriptFolder);
addpath(p1);
structFolder = [resultsFolder,'ptStructs/'];

seq_file = 'long_seq.mat';
power_file = 'power2.mat';

load([structFolder,seq_file])

if exist([structFolder,power_file],'file') ~= 0
    load([structFolder,power_file])
    
else
    power = struct;
end


if isempty(whichPts) == 1
    for i = 1:length(pt)
        if isempty(pt(i).seq_matrix) == 0
            whichPts = [whichPts,i];
        end
    end
end

else
    example = 1;
    power = struct;
end

for whichPt = whichPts
    fprintf('Doing %s\n',pt(whichPt).name)
    
    %{
    if size(power,1) >= whichPt
        if isfield(power(whichPt),'ad_rat_fft') == 1
            if isempty(power(whichPt).ad_rat_fft) == 0
                fprintf('Already did %s, skipping\n',pt(whichPt).name);
                continue
            end
        end
    else
    
    end
    %}
    if example == 0
        dataName = pt(whichPt).ieeg_name;
    else
        dataName = [];
        pwfile = [];
    end
    fs = pt(whichPt).fs;
    channels = pt(whichPt).channels;
    nch = length(channels);
    
    nbands = ceil(fs/2/10);
    
    szTimes = pt(whichPt).newSzTimes;
    
    power(whichPt).times = mean(pt(whichPt).runTimes,2);
    if isfield(power(whichPt),'alpha') == 0 || isempty(power(whichPt).alpha) == 1
        power(whichPt).alpha = zeros(nch,size(pt(whichPt).runTimes,1));
        power(whichPt).delta = zeros(nch,size(pt(whichPt).runTimes,1));
        
        %power(whichPt).ad_rat = zeros(nch,size(pt(whichPt).runTimes,1));
        %power(whichPt).ad_rat_alt = zeros(nch,size(pt(whichPt).runTimes,1));
        
        power(whichPt).finished = zeros(size(pt(whichPt).runTimes,1),1);
    end
    %  Loop over run times
    for tt = 1:size(pt(whichPt).runTimes,1)
        
        if power(whichPt).finished(tt) == 1
            fprintf('Already did chunk %d for %s, skipping\n',...
                tt,pt(whichPt).name);
            continue
        end
        
        fprintf('Doing chunk %d of %d for %s\n',tt,...
            size(pt(whichPt).runTimes,1),pt(whichPt).name);
        
        
        if example == 0
        % Add a button push to the desmond file (for the purpose of
        % restarting the program if it crashes due to java heap errors (I
        % think the ieeg toolbox creates a memory leak...))
        buttonpush = datestr(now,'yyyy-mm-dd HH:MM:SS');
        allwrite = [buttonpush,'\n'];
        fid = fopen('/tmp/desmond.txt','wt');
        fprintf(fid,allwrite);
        fclose(fid);
        end
        
        
        % Get the desired indices
        desiredTimes = pt(whichPt).runTimes(tt,:);
        
        % get times to clip
        clipTime = [-1*60 0];
        szTimesPlusClip = szTimes + repmat(clipTime,size(szTimes,1),1);
        szTimesT = szTimesPlusClip';
        out=range_intersection(desiredTimes,szTimesT(:));
        if isempty(out) == 1
            indicesToClip = []; 
        else
            indicesToClip = round(out(1)*fs):round(out(2)*fs);
        end
        
        %desiredTimes = desiredTimesFake(tt,:);
        indices = max(round(desiredTimes(1)*fs),1):round(desiredTimes(2)*fs);
        
        
        % Run the code that actually calculates the alpha delta ratio
        power(whichPt).alpha(:,tt) = innerAlphaDelta(dataName,channels,indices,pwfile,indicesToClip,fs,3,pt,whichPt);
        power(whichPt).delta(:,tt) = innerAlphaDelta(dataName,channels,indices,pwfile,indicesToClip,fs,4,pt,whichPt);
        power(whichPt).finished(tt) = 1;
        if example == 0
        save([structFolder,power_file],'power')
        end
        %fprintf('Finished analysis\n');
        
        
       
    
    end
    
  %  power(whichPt).all_p = all_p;
   % power(whichPt).ad_rat_fft = ad_rat;
   % power(whichPt).ad_rat_band = ad_rat_band;
    %power(whichPt).times = times_out;
    
end

if example == 0
% Make a new document if I make it here
fid2 = fopen('/tmp/ok.txt','wt');
fprintf(fid2,'Done\n');
%fflush(fid2);
fclose(fid2);
end


end