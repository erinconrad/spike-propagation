function icInteric(pt_ic,pt,whichPt)

% output file name
[~,~,~,resultsFolder,~] = fileLocations;

% Get wij
xyChan = pt(whichPt).electrodeData.locs;
wij = getwij(xyChan,pt(whichPt).dmin);
nchs = length(pt(whichPt).channels);

% Get interictal sequences
[seq_inter,times_inter] = divideIntoSzChunks(pt,whichPt);
seq_inter(seq_inter==0) = nan; % WHY ARE THERE ANY ZEROS?????

% Get ictal sequences
seq_ic = [];
for j = 1:length(pt_ic(whichPt).sz)
    seq_ic = [seq_ic,pt_ic(whichPt).sz(j).seq_matrix];
end
times_ic = min(seq_ic,[],1);



% get average recruitment latency for ictal and interictal sequences
%[RL,~] = getRL(seqs);
[RL_ic,~] = getRL(seq_ic);
[RL_inter,~] = getRL(seq_inter);

% Get MI
MIstruct_ic = moranStats(RL_ic',wij,nchs);
MIstruct_inter = moranStats(RL_inter',wij,nchs); 

% Get vectors (out of order!)
[vec_all,early_all,late_all] = getVectors2([seq_ic,seq_inter],pt(whichPt).electrodeData);
[vec_ic,early_ic,late_ic] = getVectors2(seq_ic,pt(whichPt).electrodeData);
[vec_inter,early_inter,late_inter] = getVectors2(seq_inter,pt(whichPt).electrodeData);

vec_ic = vec_ic./vecnorm(vec_ic,2,2);
vec_inter = vec_inter./vecnorm(vec_inter,2,2);
vec_all = vec_all./vecnorm(vec_all,2,2);

% Plot vectors over time

smspan = 1;
%{
figure
subplot(2,1,1)
scatter(times_inter,smooth(vec_inter(:,1),smspan),'b');
yl=ylim;
hold on
scatter(times_inter,smooth(vec_inter(:,2),smspan)+2,'r');
yl=ylim;
scatter(times_inter,smooth(vec_inter(:,3),smspan)+4,'g');

subplot(2,1,2)
scatter(1:length(vec_ic),smooth(vec_ic(:,1),smspan),'b');
yl=ylim;
hold on
scatter(1:length(vec_ic),smooth(vec_ic(:,2),smspan)+2,'r');
yl=ylim;
scatter(1:length(vec_ic),smooth(vec_ic(:,3),smspan)+4,'g');


figure
subplot(2,1,1)
scatter(times_inter,smooth(early_inter(:,1),smspan),'b');
yl=ylim;
hold on
scatter(times_inter,smooth(early_inter(:,2),smspan)+2,'r');
yl=ylim;
scatter(times_inter,smooth(early_inter(:,3),smspan)+4,'g');

subplot(2,1,2)
scatter(1:length(vec_ic),smooth(early_ic(:,1),smspan),'b');
yl=ylim;
hold on
scatter(1:length(vec_ic),smooth(early_ic(:,2),smspan)+2,'r');
yl=ylim;
scatter(1:length(vec_ic),smooth(early_ic(:,3),smspan)+4,'g');
%}

all_vecs= [vec_ic;vec_inter];
icIdx = [ones(size(vec_ic,1),1);zeros(size(vec_inter,1),1)];

% MANOVA comparing vectors for ictal and interictal sequences
[d,p,stats] = manova1(all_vecs,icIdx,0.01);

% Get first sequences channels
[~,firstChsInter] = min(seq_inter,[],1);
firstChLocsInter = xyChan(firstChsInter,2:4);
[~,firstChsIc] = min(seq_ic,[],1);
firstChLocsIc = xyChan(firstChsIc,2:4);

firstLocsAll = [firstChLocsIc;firstChLocsInter];
icIdx = [ones(size(firstChLocsIc,1),1);zeros(size(firstChLocsInter,1),1)];

% clustering algorithm
[idx,C,sumd,D] = kmeans([firstLocsAll,vec_all],3);
colors = [1 0 0;0 1 0;0 0 1; 0.5 0.5 1; 1 0.5 0.5; 0.5 1 0.5];

c_idx = zeros(size(idx,1),3);
for i = 1:length(idx)
   c_idx(i,:) = colors(idx(i),:); 
    
end

figure
subplot(3,1,1)
scatter(1:length(vec_ic),smooth(firstChLocsIc(:,1),smspan),20,c_idx(logical(icIdx),:),'o');
yl=ylim;
hold on
scatter(1:length(vec_ic),smooth(firstChLocsIc(:,2)+100,smspan),20,c_idx(logical(icIdx),:),'x');
yl=ylim;
scatter(1:length(vec_ic),smooth(firstChLocsIc(:,3),smspan),20,c_idx(logical(icIdx),:),'*');

subplot(3,1,2)
scatter(1:length(vec_ic),smooth(idx(logical(icIdx)),smspan),20,c_idx(logical(icIdx),:));

subplot(3,1,3)
scatter3(xyChan(:,2),xyChan(:,3),xyChan(:,4),60,'k');
hold on
for k = 1:size(C,1)
    scatter3(C(k,1),C(k,2),C(k,3),60,colors(k,:),'filled');
    plot3([C(k,1) C(k,1) + C(k,4)],...
        [C(k,2) C(k,2) + C(k,5)],...
        [C(k,3) C(k,3) + C(k,6)],'k','LineWidth',2)
end

figure
subplot(3,1,1)
scatter(times_inter/3600,smooth(firstChLocsInter(:,1),smspan),20,c_idx(~logical(icIdx),:),'o');
hold on
scatter(times_inter/3600,smooth(firstChLocsInter(:,2)+100,smspan),20,c_idx(~logical(icIdx),:),'x');
scatter(times_inter/3600,smooth(firstChLocsInter(:,3),smspan),20,c_idx(~logical(icIdx),:),'*');
for j = 1:length(pt(whichPt).sz)
   yl = ylim; 
   szOnset = pt(whichPt).sz(j).onset;
   sz = plot([szOnset szOnset]/3600,yl,'k--','LineWidth',2);
end

subplot(3,1,2)
scatter(times_inter/3600,smooth(idx(~logical(icIdx)),smspan),20,c_idx(~logical(icIdx),:));

subplot(3,1,3)
scatter3(xyChan(:,2),xyChan(:,3),xyChan(:,4),60,'k');
hold on
for k = 1:size(C,1)
    scatter3(C(k,1),C(k,2),C(k,3),60,colors(k,:),'filled');
    plot3([C(k,1) C(k,1) + C(k,4)],...
        [C(k,2) C(k,2) + C(k,5)],...
        [C(k,3) C(k,3) + C(k,6)],'k','LineWidth',2)
end



%{

figure
subplot(2,1,1)
scatter(times_inter,smooth(firstChLocsInter(:,1),smspan),'b');
yl=ylim;
hold on
scatter(times_inter,smooth(firstChLocsInter(:,2),smspan)+2,'r');
yl=ylim;
scatter(times_inter,smooth(firstChLocsInter(:,3),smspan)+4,'g');

subplot(2,1,2)
scatter(1:length(vec_ic),smooth(firstChLocsIc(:,1),smspan),'b');
yl=ylim;
hold on
scatter(1:length(vec_ic),smooth(firstChLocsIc(:,2),smspan)+2,'r');
yl=ylim;
scatter(1:length(vec_ic),smooth(firstChLocsIc(:,3),smspan)+4,'g');
%}

% the next step will be to figure out how to combine it for multiple
% seizures and multiple patients with manova
% https://www.mathworks.com/help/stats/repeatedmeasuresmodel.manova.html


% Test that the vectors form more or less a normal distribution
%{
mean_vec = mean(vec_all,1);
std_vec = std(vec_all,1);
n_vecs = size(vec_all,1);
z = repmat(mean_vec,n_vecs,1) + randn(n_vecs,3).*std_vec;
scatter3(z(:,1),z(:,2),z(:,3),'g');
hold on
scatter3(vec_all(:,1),vec_all(:,2),vec_all(:,3),'b');
%}

% Compare the ictal and interictal vectors spatially
%{
figure
scatter3(vec_ic(:,1),vec_ic(:,2),vec_ic(:,3),'b')
hold on
scatter3(vec_inter(:,1),vec_inter(:,2),vec_inter(:,3),'r')

%}



end