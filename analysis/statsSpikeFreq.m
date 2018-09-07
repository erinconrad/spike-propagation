function statsSpikeFreq(pt,whichPts)

%% Stats plan
% First I will try an n-way ANOVA, where the dependent variable is spike
% frequency and the independent variables are 1) which patient, 2) which
% seizure, and 3) which time chunk preceding the seizure

% I will only look at seizures that are at least 6 hours after the last
% seizure and have at least 6 hours of spike data preceding them 

% Consider removing the last 30 minutes or so before the seizure


%% Remove EKG artifact and depth electrodes
rmEKG = 1;
prox = 0.02; %20 ms
rmDepth = 1;
rmType = 'D';

%% Parameters
window = 3600;
nwindows = 6;


[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
gdfFolder = [resultsFolder,'gdf/'];


%% Initialize stats variables
whichPt = [];
whichSz = [];
whichChunk = [];
spikeNum = [];
szNum = 0;


for i = whichPts
    for j = 1:length(pt(i).sz)
        
        gdf_all = [];
        
        %% Get all the spikes for the seizure
        
        % Skip if too close to the beginning of the file
        if pt(i).sz(j).onset < window*nwindows
            continue
        end
        
        % Skip if too close to the last seizure
        if j > 2 && pt(i).sz(j).onset - pt(i).sz(j-1).onset < window*nwindows
            continue
        end
        
        for k = 1:length(pt(i).sz(j).chunkFiles)
            
            if exist([gdfFolder,pt(i).name,'/',pt(i).sz(j).chunkFiles{k}],'file') == 0
                continue
            end

            % Load gdf file
            load([gdfFolder,pt(i).name,'/',pt(i).sz(j).chunkFiles{k}]);

            if isempty(gdf) == 1
                continue
            end
            
            % Skip it if the last spike is more than 6 hours before the
            % seizure
            if pt(i).sz(j).onset - gdf(end,2) > window*nwindows
                continue
            end
            
            % Break the loop if the first spike is after the seizure
            if gdf(1,2) >= pt(i).sz(j).onset
                break
            end
            
            % Load gdf ekg file
            load([gdfFolder,pt(i).name,'/',pt(i).sz(j).EKGchunkFiles{k}]);

            % remove EKG artifact
            if rmEKG == 1
                gdf = removeEKGArtifact(gdf,gdf_ekg,prox);
            end

            if rmDepth == 1
                gdf = removeChs(gdf,pt(i).electrodeData,rmType);
            end
            
            gdf_all = [gdf_all;gdf];
            
        end
        
        %% Divide the spikes into the appropriate chunks
        
        % Skip this seizure if the first spike is within 5 hours of the
        % seizure (happens when there are periods of disconnection, e.g.
        % HUP070)
        if pt(i).sz(j).onset - gdf_all(1,2) < window*(nwindows-1)
            continue
        end
        
        % Increment the seizure number
        szNum = szNum + 1;
        
        % Remove spikes more than 6 hours
        gdf_all(pt(i).sz(j).onset - gdf_all(:,2) > window*nwindows ,:) = [];
        
        % Remove spikes occuring after the seizure onset
        gdf_all(gdf_all(:,2) > pt(i).sz(j).onset,:) = [];
        
        % Loop through the windows
        for tt = 1:nwindows
           times = [pt(i).sz(j).onset - window*nwindows + (tt-1)*window,...
               pt(i).sz(j).onset - window*nwindows + (tt)*window];
           
           % Get the appropriate spikes in that window
           gdf_chunk = gdf_all(gdf_all(:,2) >= times(1) & gdf_all(:,2) <= ...
               times(2),:);
           
           % Add this data to the spike frequency array and the arrays
           % containing info on the independent variables for the anova
           spikeNum = [spikeNum,size(gdf_chunk,1)];
           whichPt = [whichPt,i];
           whichSz = [whichSz,szNum];
           whichChunk = [whichChunk,tt];
            
        end


    end

end

%% Poisson regression
[b,dev,stats] = glmfit([whichChunk;whichSz]',spikeNum,'poisson','link','log');

%% Run the ANOVA
%{
[p,tbl,stats] = anovan(spikeFreq,{whichPt, whichSz, whichChunk},'model',...
    'interaction','varnames',{'whichPt','whichSz','whichChunk'});

[p,tbl,stats] = anovan(spikeFreq,{whichSz, whichChunk},'model',...
    'interaction','varnames',{'whichSz','whichChunk'});

[p,tbl,stats] = anova(spikeFreq,whichChunk)
%}

end