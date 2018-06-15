function [allSens,allAcc]=loopSpikeChecker(whichDetector,trainOrTest)

if trainOrTest == 2
    error('Are you sure you want to look at the testing data?\n');
end

tmuls_to_try = [11,12,13,14];
absthresh_to_try = [300];

%% File names
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
%electrodeFile = [electrodeFolder,csvFile];
p1 = genpath(scriptFolder);
addpath(p1);
timeFile = 'ptWithElectrodeData.mat'; 
newptfile = 'ptAccuracies.mat';
validatedFile = 'validated.mat';

%% Load file with filenames and run times
load([resultsFolder,'ptStructs/',timeFile]);
load([resultsFolder,'validation/',validatedFile]);

for i = 1:length(validated)
    
    if isempty(pt(i).ieeg_name) == 1
        fprintf('Missing ieeg_name for patient %s, skipping\n',pt(i).name);
        continue
    end
   
    if strcmp(pt(i).name,validated(i).name)~=1
        error(sprintf('Warning, pt struct ids do not match up for patient %d\n',i));
    end
    
    mkdir([resultsFolder,'validation/',pt(i).name]);
    mkdir([resultsFolder,'validation/',pt(i).name,'/train']);
    mkdir([resultsFolder,'validation/',pt(i).name,'/test']);
    
     chnames = validated(i).chs;
     
     if trainOrTest == 1
         spikeTimes = validated(i).spike_times(validated(i).train);
         notSpikeTimes = validated(i).not_spike_times(validated(i).train);
         whichSpikes = validated(i).train;
         outputDest = [resultsFolder,'validation/',pt(i).name,'/train/sensAndAccs.mat'];
     elseif trainOrTest == 2
         spikeTimes = validated(i).spike_times(validated(i).test);
         notSpikeTimes = validated(i).not_spike_times(validated(i).test);
         whichSpikes = validated(i).test;
         outputDest = [resultsFolder,'validation/',pt(i).name,'/test/sensAndAccs.mat'];

     end
    
    if isempty(spikeTimes) || isempty(chnames)
        fprintf('Missing start times or channel names for patient %s\n',P(i).name);
        continue
    end
    
    chIds = zeros(size(chnames));
    for ich = 1:length(chnames)
        [Lia,chIds(ich)] = ismember(chnames{ich},pt(i).electrodeData.unignoredChs);
        if Lia == 0
            fprintf('Warning, could not find channel %s in the list of unignored channels for patient %s\n',...
                chnames{ich},pt(i).name);
            error('');
        end
    end
       
    allSens = [];
    allAcc = [];
    
    for k = tmuls_to_try

        for m = absthresh_to_try

            [sensitivity,accuracy] = spikeChecker(pt,i,chIds,...
   spikeTimes,notSpikeTimes,k,m,whichDetector,trainOrTest,whichSpikes);
            
            
           
            
            allSens = [allSens;k m sensitivity];
            allAcc = [allAcc;k m accuracy];

        end



    end
        
    
    
%% Save output file
save(outputDest,'allSens','allAcc');
    
    
    
end


end


