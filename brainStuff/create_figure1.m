%{

- reconstruction freesurfer; renderings blender
- co-registration: 
    - electrode coordinates determined in CT, coregistered
into T1 space 
    - align electrode coordinates with MRI in blender
    - 

%}


% HUP68 offset = [-1.0029,2.3087,28.9465];

close all
clc

addpath(genpath('~/Documents/gdrive/NODDI/util/NIfTI_20140122'))
addpath(genpath('~/Documents/gdrive/NODDI/util/gifti-1.6'))

DATA_DIR = '/gdrive/public/USERS/lkini/aim1/data';
FIG_DIR = '/gdrive/public/USERS/erinconr/projects/spike_propagation/scripts/spike-propagation/results/brainPlots/';
COMP_DIR = '/gdrive/public/USERS/lkini/aim1/results';
%PTX = {'HUP052','HUP054','HUP065','HUP068','HUP078','HUP082','HUP086','HUP093','3T_P16','7T_P12'};
% PTX = {'HUP093'};
PTX = {'HUP078'};
for i = 1:length(PTX)
    PATIENT_ID = PTX{i};    
    RESULTS_FOLDER = fullfile(COMP_DIR,PATIENT_ID,'surf');
    
    iter = 1;
    filePattern = fullfile(RESULTS_FOLDER,'lh*FLAIR_2*JACOBIAN-in-brain.gii');
    files = dir(filePattern);
    resection_filePattern = fullfile(RESULTS_FOLDER,'lh*RESECTION_final_with_mask.gii');
    resection_file = dir(resection_filePattern);
    g = gifti(fullfile(RESULTS_FOLDER,'lh.pial.gii'));
    for ii = 1:length(files)
        baseFileName = files(ii).name;
        fullFileName = fullfile(RESULTS_FOLDER,baseFileName);        
        gg = gifti(fullFileName);        
        resection_fullFileName = fullfile(RESULTS_FOLDER,resection_file(1).name);
        resection_gg = gifti(resection_fullFileName);
        gg.cdata(resection_gg.cdata == 0) = zscore(gg.cdata(resection_gg.cdata == 0));
        gg.cdata(resection_gg.cdata > 0) = -100;
        h = figure(iter); plot(g,gg), caxis([-5 5]), view(-90, 0);
        saveas(h,fullfile(FIG_DIR,sprintf('%s_%i.png',PATIENT_ID,iter)));
        close all;
        iter = iter+1;
    end
    
    iter = 1;
    filePattern = fullfile(RESULTS_FOLDER,'rh*FLAIR_2*JACOBIAN-in-brain.gii');
    files = dir(filePattern);
    resection_filePattern = fullfile(RESULTS_FOLDER,'rh*RESECTION_final_with_mask.gii');
    resection_file = dir(resection_filePattern);
    g = gifti(fullfile(RESULTS_FOLDER,'rh.pial.gii'));
    for ii = 1:length(files)
        baseFileName = files(ii).name;
        fullFileName = fullfile(RESULTS_FOLDER,baseFileName);        
        gg = gifti(fullFileName);        
        resection_fullFileName = fullfile(RESULTS_FOLDER,resection_file(1).name);
        resection_gg = gifti(resection_fullFileName);
        gg.cdata(resection_gg.cdata == 0) = zscore(gg.cdata(resection_gg.cdata == 0));
        gg.cdata(resection_gg.cdata > 0) = -100;
        h = figure(iter); plot(g,gg), caxis([-5 5]), view(90, 0);
        saveas(h,fullfile(FIG_DIR,sprintf('%s_%i.png',PATIENT_ID,iter+100)));
        close all;
        iter = iter+1;
    end
    
    iter = 1;
    filePattern = fullfile(RESULTS_FOLDER,'lh*FLAIR_2*JACOBIAN-in-brain.gii');
    files = dir(filePattern);
    resection_filePattern = fullfile(RESULTS_FOLDER,'lh*RESECTION_final_with_mask.gii');
    resection_file = dir(resection_filePattern);
    g = gifti(fullfile(RESULTS_FOLDER,'lh.pial.gii'));
    for ii = 1:length(files)
        baseFileName = files(ii).name;
        fullFileName = fullfile(RESULTS_FOLDER,baseFileName);        
        gg = gifti(fullFileName);
        gg.cdata = -gg.cdata;        
        resection_fullFileName = fullfile(RESULTS_FOLDER,resection_file(1).name);
        resection_gg = gifti(resection_fullFileName);
        gg.cdata(resection_gg.cdata == 0) = zscore(gg.cdata(resection_gg.cdata == 0));
        gg.cdata(resection_gg.cdata > 0) = -100;
        h = figure(iter); plot(g,gg), caxis([-5 5]), view(-90, 0);
        saveas(h,fullfile(FIG_DIR,sprintf('n%s_%i.png',PATIENT_ID,iter)));
        close all;
        iter = iter+1;
    end
    
    iter = 1;
    filePattern = fullfile(RESULTS_FOLDER,'rh*FLAIR_2*JACOBIAN-in-brain.gii');
    files = dir(filePattern);
    resection_filePattern = fullfile(RESULTS_FOLDER,'rh*RESECTION_final_with_mask.gii');
    resection_file = dir(resection_filePattern);
    g = gifti(fullfile(RESULTS_FOLDER,'rh.pial.gii'));
    for ii = 1:length(files)
        baseFileName = files(ii).name;
        fullFileName = fullfile(RESULTS_FOLDER,baseFileName);        
        gg = gifti(fullFileName);
        gg.cdata = -gg.cdata;        
        resection_fullFileName = fullfile(RESULTS_FOLDER,resection_file(1).name);
        resection_gg = gifti(resection_fullFileName);
        gg.cdata(resection_gg.cdata == 0) = zscore(gg.cdata(resection_gg.cdata == 0));
        gg.cdata(resection_gg.cdata > 0) = -100;
        h = figure(iter); plot(g,gg), caxis([-5 5]), view(90, 0);
        saveas(h,fullfile(FIG_DIR,sprintf('n%s_%i.png',PATIENT_ID,iter+100)));
        close all;
        iter = iter+1;
    end
end