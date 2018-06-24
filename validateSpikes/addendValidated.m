function addendValidated(whichPt,whichPtName)

%{
This function allows me to add additional validated true positive and true
negative spikes to the "validated.mat" file that I will reference when
checking the accuracy of the spike detector
%}

%% File names
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
%electrodeFile = [electrodeFolder,csvFile];
p1 = genpath(scriptFolder);
addpath(p1);
validatedFile = 'validated.mat';

%% Load file with filenames and run times
a = load([resultsFolder,'validation/',validatedFile]);
validated = a.validated;

validated(whichPt) = trueSpikes(whichPtName);


%% Save the file
save([resultsFolder,'validation/',validatedFile],'validated');

end