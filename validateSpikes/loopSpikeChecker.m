function loopSpikeChecker(whichPts,whichDetector,trainOrTest,merge,tmuls_to_try,absthresh_to_try)

if trainOrTest == 2
    error('Are you sure you want to look at the testing data?\n');
end

% If merge is 1, then I will not overwrite prior checks, just add new ones
% to the overall sensitivity and accuracy. If it is 0, I will ignore prior
% checks and overwrite everything.

%% File names
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
%electrodeFile = [electrodeFolder,csvFile];
p1 = genpath(scriptFolder);
addpath(p1);
timeFile = 'ptWithElectrodeData.mat'; 
newptfile = 'ptAccuracies.mat';
validatedFile = 'validated.mat';

%% Load file with filenames and run times
load([resultsFolder,'ptStructs/allPtsElectrodeData/',timeFile]);
load([resultsFolder,'validation/',validatedFile]);



for i = whichPts
    
    if isempty(validated(i).train) == 1
        fprintf('Missing training data for patient %s, skipping \n',pt(i).name);
        continue
    end
    
    if isempty(pt(i).ieeg_name) == 1
        [pt(i).ieeg_name,pt(i).electrode_name,thresh] =  ieegAndElectodeNames(pt(i).name);
        if isempty(pt(i).ieeg_name) == 1
            fprintf('Missing ieeg_name for patient %s, skipping\n',pt(i).name);
            continue
        end
    end
    
    if isempty(tmuls_to_try) == 1
        tmuls_to_try =  thresh.tmul;
    end
    
    if isempty(absthresh_to_try) == 1
        absthresh_to_try = thresh.absthresh;
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
         outputDest = [resultsFolder,'validation/',pt(i).name,'/train/sensAndAccsDet',sprintf('%d',whichDetector),'.mat'];
     elseif trainOrTest == 2
         spikeTimes = validated(i).spike_times(validated(i).test);
         notSpikeTimes = validated(i).not_spike_times(validated(i).test);
         whichSpikes = validated(i).test;
         outputDest = [resultsFolder,'validation/',pt(i).name,'/test/sensAndAccsDet',sprintf('%d',whichDetector),'.mat'];

     end
     
     
    
    
    if isempty(spikeTimes) || isempty(chnames)
        fprintf('Missing start times or channel names for patient %s\n',pt(i).name);
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
       
    
    
    for k = tmuls_to_try

        for m = absthresh_to_try
            
            
            
            fprintf('Doing tmul %d and absthresh %d for patient %s\n',k,m,pt(i).name);
            
             %% Try to load struct with accuracies if it exists
            if merge == 1 && exist(outputDest,'file') ~= 0
                a = load(outputDest);
                oldAllSens = a.allSens; oldAllAcc = a.allAcc;
            end
            
            if merge == 1 && exist(outputDest,'file') ~= 0
               if sum(ismember([k,m],oldAllSens)) == 2
                  fprintf('Already did tmul %d and absthresh %d, skipping...\n',k,m);
                  continue; 
                   
               end
                
            end
            
            % Add a button push to the desmond file (for the purpose of
            % restarting the program if it crashes due to random server
            % error)
            buttonpush = datestr(now,'yyyy-mm-dd HH:MM:SS');
            allwrite = [buttonpush,'\n',sprintf('Patient %s tmul %d absthresh %d\n',...
                i,k,m)];
            fid = fopen('/tmp/desmond_valid.txt','wt');
            fprintf(fid,allwrite);
            fclose(fid);

            [sensitivity,accuracy] = spikeChecker(pt,i,chIds,...
   spikeTimes,notSpikeTimes,k,m,whichDetector,trainOrTest,whichSpikes);
            
            
           
            
            allSens = [k m sensitivity];
            allAcc = [k m accuracy];
            
            %% Save output file
            if merge == 1 && exist(outputDest,'file') ~= 0
                allSens = unique([oldAllSens;allSens],'rows');
                allAcc = unique([oldAllAcc;allAcc],'rows');
                
                
            end
            
            if merge == 1
                save(outputDest,'allSens','allAcc');

                validated(i).allSens = allSens;
                validated(i).allAcc = allAcc;

                save([resultsFolder,'validation/validated.mat'],'validated');
            end
            
            

        end



    end
        
    
    

    
    
end




end


