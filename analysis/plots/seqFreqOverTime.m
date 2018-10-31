function [chunk_seqs,times_plot,MI,rl,angle,chunk_seqs_chs,vec,early] = seqFreqOverTime(pt,whichPt,window)


%% Parameters
dmin = pt(whichPt).dmin;
nchs = length(pt(whichPt).channels);

% output file name
[~,~,~,resultsFolder,~] = fileLocations;

%% Get wij
xyChan = pt(whichPt).electrodeData.locs;
wij = getwij(xyChan,dmin);

%% First need to divide it up into multiple seizure chunks for plotting

[all_seq_cat,all_times,seq_all,trackingNo] = divideIntoSzChunksGen(pt,whichPt);


chunk_seqs = cell(size(seq_all));
chunk_seqs_chs = cell(size(seq_all));
times_plot = cell(size(seq_all));
rl = cell(size(seq_all));
MI = cell(size(seq_all));
vec = cell(size(seq_all));
dot_prod = cell(size(seq_all));
angle = cell(size(seq_all));

%% Get a reference vector

%firstSpikes = min(seq,[],1);
%first_sz_seq = seq(:,firstSpikes >= pt(whichPt).sz(1).onset...
%    & firstSpikes <= pt(whichPt).sz(1).offset);
ref_vec = mean(getVectors2(all_seq_cat,pt(whichPt).electrodeData));


%% REMOVE ME
%{
ref_vec = pt(whichPt).electrodeData.ref_vector(2,:)-...
pt(whichPt).electrodeData.ref_vector(1,:);
%}

all_angle = [];

%% Now divide the sequences into windows
for i = 1:length(seq_all)
   seq = seq_all{i};
   firstSpikes = min(seq,[],1);
   totalTime = firstSpikes(end) - firstSpikes(1);
   nchunks = ceil(totalTime/window);
   
   chunk_seqs{i} = zeros(nchunks,1);
   chunk_seqs_chs{i} = zeros(nchunks,nchs);
   times_plot{i} = zeros(nchunks,1);
   rl{i} = zeros(nchunks,nchs);
   MI{i} = zeros(nchunks,1);
   vec{i} = zeros(nchunks,3);
   dot_prod{i} = zeros(nchunks,1);
   early{i} = zeros(nchunks,3);
   
   for tt = 1:nchunks
      times =  [(tt-1)*window + firstSpikes(1),tt*window + firstSpikes(1)];
      times_plot{i}(tt) = (times(1)+times(2))/2;
      
      % Get the appropriate sequences in this time
      correct_seqs = seq(:,firstSpikes >= times(1) & firstSpikes <= times(2));
    
      chunk_seqs{i}(tt) = size(correct_seqs,2);
      
      % get the number of sequences per channel (the starting channel???)
      for k = 1:size(correct_seqs,2)
         [~,ch] =  min(correct_seqs(:,k));
         chunk_seqs_chs{i}(tt,ch) = chunk_seqs_chs{i}(tt,ch) + 1;
          
      end
      
      %% Get recruitment latency for these sequences
      
      % for each sequence, the latency with which each channel is activated
      % in the sequence is the spike time in that channel minus the spike
      % time in the channel activated the earliest in that sequence
      latency_all_seq = correct_seqs - min(correct_seqs,[],1);
      
      % Take the average latency for the channel over all sequences
      mean_latency = nanmean(latency_all_seq,2);
      rl{i}(tt,:) = mean_latency';
      
      % Get the moran index
      MIstruct= moranStats(mean_latency',wij,nchs);
      if MIstruct.I > 1
          %error('look\n');
      end
      MI{i}(tt) = MIstruct.I;
      
      %% Get a vector representing the direction of each sequence
      % Confirm that this function is correct!!!!!
      [vec_temp,early_temp] = getVectors2(correct_seqs,pt(whichPt).electrodeData);
      vec{i}(tt,:) = mean(vec_temp,1);
      early{i}(tt,:) = mean(early_temp,1);
      dot_prod{i}(tt) = dot(mean(vec_temp,1)/norm(mean(vec_temp,1)),...
          ref_vec/norm(ref_vec));
      angle{i}(tt) = acos(dot_prod{i}(tt))*180/pi;

      % Get individual angle
%{
      dots = dot(vec_temp/norm(vec_temp),...
    repmat(ref_vec/norm(ref_vec),size(vec_temp,1),1),2);
      angles = acos(dots)*180/pi;
      all_angle = [all_angle;angles];
%}
      
   end
   
   chunk_seqs{i} = chunk_seqs{i}/window;
   chunk_seqs_chs{i} = chunk_seqs_chs{i}/window;
   
    
end


%% Look at all angles
%{

histogram(all_angle);

fake = normrnd(mean(all_angle),std(all_angle),size(all_angle));

figure
histogram(fake);
%}

%{
figure
for i = 1:length(chunk_seqs)
    times = times_plot{i};
    plot(times/3600,chunk_seqs{i},'k');
    hold on
end

yl = ylim;

for j = 1:length(pt(whichPt).sz)
    szTimes = [pt(whichPt).sz(j).onset,pt(whichPt).sz(j).offset];
    meanSzTimes = (szTimes(1) + szTimes(2))/2;
    plot([meanSzTimes meanSzTimes]/3600,[yl(1) yl(2)],'k--');
end

xlabel('Hour');
ylabel('Sequences per hour');
title('Sequence frequency over time');
set(gca,'FontSize',15)
%}


end