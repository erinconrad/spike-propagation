function glme = statsSpikeFreq(pt,whichPts,window)

%% Remaining problems
% I have no way to remove seizures from the analysis if there are large
% chunks of time in the 6 hours prior to the seizure where data is missing


%% Remove EKG artifact and depth electrodes
rmEKG = 1;
prox = 0.02; %20 ms
rmDepth = 1;
rmType = 'D';

%% Parameters
nwindows = 6*3600/window;

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
gdfFolder = [resultsFolder,'gdf/'];


%% Initialize stats variables

% A cell that is of size 1xn, where n is the number of retained seizures.
% Each element in the cell is the gdf_all for that seizure
gdf_diff_sz = {};

% Array of all spike times
spikes = [];

% The number of the seizure index
szIndex = 0;

% An mx1 array of all spike counts where m is the number of windows
all_num_spikes = [];

% An mx1 array of all patient numbers for each window where m is the number
% of windows
pat = [];

% mx1 array of all seizure index numbers
sz = [];

% mx1 array of all window numbers
chunk = [];


% Loop through the patients
for i = whichPts
    
    % Loop through all the seizures for each patient
    for j = 1:length(pt(i).sz)
        
        % initialize the gdf for that seizure
        gdf_all = [];
        
        %% Get all the spikes for the seizure
        
        % Skip if too close to the beginning of the file
        if pt(i).sz(j).onset < window*nwindows
            continue
        end
        
        % Skip if too close to the previous seizure
        if j > 2 && pt(i).sz(j).onset - pt(i).sz(j-1).onset < window*nwindows
            continue
        end
        
        
        
        % Loop through all the gdf files for the seizure
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
            
            % Add the gdf to the full gdf for the seizure
            gdf_all = [gdf_all;gdf];
            
        end
        
        
        % Skip this seizure if the first spike is within 5 hours of the
        % seizure (happens when there are periods of disconnection, e.g.
        % HUP070)
        if pt(i).sz(j).onset - gdf_all(1,2) < window*(nwindows-1)
            continue
        end
        
        % We are going to use this seizure, so increase the seizure index
        % by 1
        szIndex = szIndex + 1;
        
        % Remove spikes more than 6 hours
        gdf_all(pt(i).sz(j).onset - gdf_all(:,2) > window*nwindows ,:) = [];
        
        % Remove spikes occuring after the seizure onset
        gdf_all(gdf_all(:,2) > pt(i).sz(j).onset,:) = [];
        
        % Add the gdf_all for this seizure to the cell array
        gdf_diff_sz{end+1} = gdf_all;
        
        % Get spike times and add it to the array of all spike times
        temp_spikes = gdf_all(:,2);
        spikes = [spikes;temp_spikes];
        
        %% Now get window data
        startTime = pt(i).sz(j).onset - window*nwindows;
        new_chunks = zeros(size(gdf_all(:,2)));
        
        % Loop through all the windows
        for tt = 1:nwindows
            
            % Get the times for the windows
            times = [(tt-1)*window + startTime,tt*window + startTime];
            
            
            new_chunks(temp_spikes>=times(1) & temp_spikes <=times(2)) = tt;
            
            % Get the number of spikes in the window
            all_num_spikes = [all_num_spikes;...
                sum(temp_spikes>=times(1) & temp_spikes <=times(2))];
            
            % Get the patient id for that window
            pat = [pat;i];
            
            % get the seizure id for that window
            sz = [sz;szIndex];
            
            % get the chunk (or window) id for that window
            chunk = [chunk;tt];
            
        end
        
        if sum(new_chunks==0) > 0
            error('look\n');
        end
        
        %{
        chunk = [chunk;new_chunks];
        pat = [pat;ones(size(new_chunks))*i];
        sz = [sz;ones(size(new_chunks))*szIndex];
        %}
        
    end

end




%% Poisson regression to combine data for all spikes and all seizures
% http://math.bu.edu/people/mak/samsi/SAMSI_GLM_Example.pdf

% This method assumes that spikes from all
% patients and all seizures come from a population with the same mean. This
% is incorrect and will likely cause problems as I include more patients
% and more seizures. However, for small numbers of seizures it appears to
% produce the same result as a generalized linear mixed model.
[b,dev,stats] = glmfit([chunk],all_num_spikes,'poisson');

%% Plot the data against the model
%{
lambda = exp(b(1) + b(2)*chunk);
figure
scatter(chunk,all_num_spikes,'b')
hold on
scatter(chunk,lambda,'r')
xlabel('which window')
ylabel('spike count')
legend('Data','GLM');
%}



%% Make a table for the purpose of doing a generalized linear mixed model
T = table(all_num_spikes,chunk,sz,pat,'VariableNames',{'all_num_spikes','chunk','sz','pat'});
glme = fitglme(T,'all_num_spikes ~ 1 + chunk + (1|sz) + (1|pat)',...
    'Distribution','Poisson','Link','log','FitMethod','Laplace',...
    'DummyVarCoding','effects');

%% Old statistics
%{

% test for each seizure if it is a homogeneous poisson process
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
    p_L = 1-normcdf(z);
    
    
    
    %% Kolmogorov-Smirnov test
    f_data = linspace(min(spike_times),max(spike_times),...
        length(spike_times));
    %f_data = f_data + 5000*rand(1,length(f_data));
    
    %{
    figure
    plot(f_data)
    hold on
    plot(spike_times,'r')
    %}
    
    fake_uniform  = makedist('uniform','lower',...
        min(spike_times),'upper',max(spike_times));
    [h,p_KS,ksstat,cv] = kstest(spike_times,'cdf',fake_uniform);

     
    %% Military handbook test
    %{
    sum_logs = 0;
    for tt = 1:length(spike_times)
    	sum_logs = sum_logs +log(spike_times(end)/spike_times(tt));
    end
    chi2 = 2*sum_logs;
    %}
    
    %http://www.kybernetika.cz/content/2007/4/415/paper.pdf
    chi2 = 2*sum(log(spike_times(end)./spike_times));
    p_MH = chi2cdf(chi2,2*n);
    
    beta_ml = n/sum(log(spike_times(end)./spike_times));
    alpha_ml = n/(spike_times(end))^beta_ml;
    %{
    figure
    plot(linspace(0,spike_times(end),n),alpha_ml*(linspace(0,spike_times(end),n)).^beta_ml,'r');
    hold on
    plot(spike_times,cum_spike_num,'b');
    plot(linspace(0,spike_times(end),n),n/spike_times(end)*linspace(0,spike_times(end),n),'g')
    title('Power law model')
    %}
    
    
    % Example of tornadoes
    %http://www.kybernetika.cz/content/2007/4/415/paper.pdf
    %{
    tornado_days = [215
        1210
        1552
        1721
        2192
        3532
        3698
        4976
        5098
        6761
        6829
        7197
        7229
        7235
        7329
        7338
        7361
        7899
        7939
        7962
        8271
        8289
        8469
        8646
        8699
        8782
        9147
        9151
        9256
        9407];
    tornado_days = tornado_days - tornado_days(1)+5;
    beta_ml_tor = length(tornado_days)/...
        sum(log(tornado_days(end)./tornado_days));
    alpha_ml_tor = length(tornado_days)/...
        (tornado_days(end))^beta_ml_tor;
    figure
    plot(tornado_days,cumsum(ones(length(tornado_days),1)),'r');
    hold on
    plot((linspace(tornado_days(1),tornado_days(end),...
        length(tornado_days))),alpha_ml_tor*(linspace(tornado_days(1),tornado_days(end),...
        length(tornado_days))).^beta_ml_tor);
    plot((linspace(tornado_days(1),tornado_days(end),...
        length(tornado_days))),...
        length(tornado_days)/tornado_days(end)*(linspace(tornado_days(1),tornado_days(end),...
        length(tornado_days))),'g');
    %}
    
end    
    
% Attempt at ANOVA
[p,tbl,stats] = anovan(spikeFreq,{whichPt, whichSz, whichChunk},'model',...
    'interaction','varnames',{'whichPt','whichSz','whichChunk'});

[p,tbl,stats] = anovan(spikeFreq,{whichSz, whichChunk},'model',...
    'interaction','varnames',{'whichSz','whichChunk'});

[p,tbl,stats] = anova(spikeFreq,whichChunk)

%}

end