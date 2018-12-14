function [p,info] = correlateClinically(measure,clinical,measureType,clinicalType,doPlot)

if strcmp(measureType,'bin') == 1 && strcmp(clinicalType,'num') == 1
    % compare ranked clinical data between two groups. Wilcoxon rank rum
    % test. (e.g., do patients who have a significant change in cluster
    % distribution between pre-ictal and inter-ictal period have different
    % clinical outcomes?)
    [p,~,stats] = ranksum(clinical(measure==1),clinical(measure==0));
    info.stats = stats;
    if doPlot == 1
        figure
        errorbar(0,mean(clinical(measure==0)),std(clinical(measure==0)),'o');
        hold on
        errorbar(1,mean(clinical(measure==1)),std(clinical(measure==1)),'o');

    end

elseif strcmp(measureType,'num') == 1 && strcmp(clinicalType,'num') == 1
    % correlate ranked clinical data and ranked measure data. Spearman rank
    % correlation. (e.g., does number of hours needed to capture cluster
    % variation correlate with outcome?)
    [rho,p] = corr(measure,clinical,'Type','Spearman');
    info.rho = rho;
    if doPlot == 1
        figure
        scatter(measure,clinical,100,'filled');
    end
elseif strcmp(measureType,'num') == 1 && strcmp(clinicalType,'bin') == 1
    % Same as the first type, wilcoxon rank sum test. (e.g., does number of
    % hours needed to capture cluster variation differ between different
    % epilepsy locations if only two possible groups?)
    [p,~,stats] = ranksum(measure(clinical==1),measure(clinical==0));
    info.stats = stats;
    if doPlot == 1
        figure
        errorbar(0,mean(measure(clinical==0)),std(measure(clinical==0)),'o');
        hold on
        errorbar(1,mean(measure(clinical==1)),std(measure(clinical==1)),'o');

    end
elseif strcmp(measureType,'num') == 1 && strcmp(clinicalType,'cat') == 1
    % Compare ranked measure data across >2 clinical groups. Kruskal wallis
    % test. (e.g., does number of hours needed to capture cluster variation
    % differ between different epilepsy locations if >2 possible
    % locations?)
    [p,tbl,stats] = kruskalwallis(measure,clinical,'off');
    info.stats = stats;
    info.tbl = tbl;
    if doPlot == 1
        figure
        for i = 1:length(unique(clinical))
            errorbar(i,mean(measure(clinical==i)),std(measure(clinical==i)),'o');
        end
    end

end


end