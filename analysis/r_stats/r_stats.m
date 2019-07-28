function r_stats(pt,cluster)

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
p1 = genpath(scriptFolder);
addpath(p1);
destFolder = [resultsFolder,'for_r/'];
mkdir(destFolder);
whichPts = [];

for i = 1:length(pt)
    if isempty(pt(i).seq_matrix) == 0
        if size(cluster(i).bad_cluster) < cluster(i).k
            whichPts = [whichPts,i];
        end
    end
end

if isequal(whichPts,[1,4,6,7,8,9,12,14,15,16,17,18,19,20,22,24,25,27,30,31]) == 0
    error('Warning, not doing correct patients!\n');
end


    
%% Load in csv files
glarma = readtable([destFolder,'glarma.csv']);

%% Find the significant patients
sig_pats = glarma.all_names(glarma.p < 0.05/20);

%% Fisher's test
% combine p values
all_p = max(glarma.p,0.001*ones(size(glarma.p)));
X_2 = -2 * sum(log(all_p));
sum_p = 1-chi2cdf(X_2,2*length(all_p));
fprintf(['There are %d with significant correlation.\nThe combined p-value'...
    'is %1.1e.\n'],sum(glarma.p<0.05/20),sum_p);


%% T test to see if the z-values are significantly different from zero
[h,p,ci,stats] = ttest(glarma.z)


%% Are the z-values for just the significant change ones consistent?
[h,p,ci,stats] = ttest(glarma.z(glarma.p<0.05/20))

%% Make pretty table

new_table = table((glarma.all_names),glarma.b,(glarma.z),glarma.p);
new_table = [new_table(1:3,:);table(('HUP075'),nan,nan,nan);new_table(4:end,:)];
new_table = [new_table(1:7,:);table(('HUP088'),nan,nan,nan);new_table(8:end,:)];
new_table = [new_table(1:9,:);table(('HUP105'),nan,nan,nan);new_table(10:end,:)];

p_sleep_t = getPText(new_table.Var4);
new_table = table(char(new_table.Var1),char(num2str(new_table.Var2,3)),...
    char(num2str(new_table.Var3,3)),char(p_sleep_t));

end


function p_text_cell = getPText(p_array)

    p_text_cell = cell(length(p_array),1);
    for i = 1:length(p_array)
        if p_array(i) < 0.001
            p_text_cell{i} = '<0.001';
        else
            p_text_cell{i} = sprintf('%1.3f',p_array(i));
        end
    
        if p_array(i) < 0.001/length(p_array)
            p_text_cell{i} = [p_text_cell{i},'***'];
        elseif p_array(i) < 0.01/length(p_array)
            p_text_cell{i} = [p_text_cell{i},'**'];
        elseif p_array(i) < 0.05/length(p_array)
            p_text_cell{i} = [p_text_cell{i},'*'];
        end
    
    end
        
        
end
    
