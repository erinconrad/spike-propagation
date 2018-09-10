function statsSpikeFreq(pt,whichPts)

%% Stats plan

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
gdf_diff_sz = {};


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
        
        
        % Skip this seizure if the first spike is within 5 hours of the
        % seizure (happens when there are periods of disconnection, e.g.
        % HUP070)
        if pt(i).sz(j).onset - gdf_all(1,2) < window*(nwindows-1)
            continue
        end
        
        
        % Remove spikes more than 6 hours
        gdf_all(pt(i).sz(j).onset - gdf_all(:,2) > window*nwindows ,:) = [];
        
        % Remove spikes occuring after the seizure onset
        gdf_all(gdf_all(:,2) > pt(i).sz(j).onset,:) = [];
        
        gdf_diff_sz{end+1} = gdf_all;
        

    end

end


for i = 1:length(gdf_diff_sz)
    gdf = gdf_diff_sz{i};
    spike_times = gdf(:,2);
    spike_times = spike_times - spike_times(1)+0.5;
    n = length(spike_times);
    cum_spike_num = cumsum(ones(size(spike_times)));
    
    %% Graph compared to exponential
    %plot(spike_times,cum_spike_num);
    %plot(
    %hold on
    
    %% Laplace test
    % https://www.itl.nist.gov/div898/handbook/apr/section2/apr234.htm#The%20Laplace
    z = sqrt(12*n)*sum(spike_times - spike_times(end)/2)/...
        (n*spike_times(end));
    p_L = 1-normcdf(z)
    
    
    %% Kolmogorov-Smirnov test
    f_data = linspace(min(spike_times),max(spike_times),...
        length(spike_times));
    %f_data = f_data + 5000*rand(1,length(f_data));
    figure
    plot(f_data)
    hold on
    plot(spike_times,'r')
    fake_uniform  = makedist('uniform','lower',...
        min(spike_times),'upper',max(spike_times));
    [h,p_KS,ksstat,cv] = kstest(spike_times,'cdf',fake_uniform)

     
    %% Military handbook test
    %{
    sum_logs = 0;
    for tt = 1:length(spike_times)
    	sum_logs = sum_logs +log(spike_times(end)/spike_times(tt));
    end
    chi2 = 2*sum_logs;
    %}
    chi2 = 2*sum(log(spike_times(end)./spike_times));
    p_MH = chi2cdf(chi2,2*n)
end





%% Poisson regression
%[b,dev,stats] = glmfit([whichChunk;whichSz]',spikeNum,'poisson','link','log');

%% Run the ANOVA
%{
[p,tbl,stats] = anovan(spikeFreq,{whichPt, whichSz, whichChunk},'model',...
    'interaction','varnames',{'whichPt','whichSz','whichChunk'});

[p,tbl,stats] = anovan(spikeFreq,{whichSz, whichChunk},'model',...
    'interaction','varnames',{'whichSz','whichChunk'});

[p,tbl,stats] = anova(spikeFreq,whichChunk)
%}

end