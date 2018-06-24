function [ xlocs,ylocs,zlocs,chanlocs,tstamps ] = spt_filegen( allseqs,xyChan );

% This function generates the necessary variables for clustering_seqs.m.
% The required inputs are xlocs, ylocs, and timestamps matrices, which
% encode the spatial and temporal information for each sequence

xyChan              = vertcat(xyChan,[0,0,0,0]);   %only changes within function
overall             = allseqs(:,1:2:end);        %toss ticks columns
overall(overall==0) = length(xyChan);            %not sure why i did this

%x and y location of each event
xlocs    = xyChan(overall(find(overall>0)),2);
ylocs    = xyChan(overall(find(overall>0)),3);
zlocs    = xyChan(overall(find(overall>0)),4);
xlocs    = reshape(xlocs,size(overall,1),size(overall,2));
ylocs    = reshape(ylocs,size(overall,1),size(overall,2));
zlocs    = reshape(zlocs,size(overall,1),size(overall,2));

%store all time stamps
overall  = allseqs(:,2:2:end);
overall(overall==0) = length(xyChan);
overall(2:end,:)    = overall(2:end,:)-repmat(overall(1,:),size(overall,1)-1,1);
overall(overall<0)  = 0;
overall(1,:)        = 0;

%return
tstamps  = overall; 
chanlocs = allseqs(:,1:2:end);

end

