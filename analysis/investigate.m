function investigate(ptOld,ptNew)

t = [40*3600 41*3600];

oldSeqs = ptOld(19).seq_matrix;
newSeqs = ptNew(19).seq_matrix;

oldTimes = min(oldSeqs,[],1);
newTimes = min(newSeqs,[],1);

%% Get sequences in this time
old_int = oldSeqs(:,oldTimes > t(1) & oldTimes < t(2));
new_int = newSeqs(:,newTimes > t(1) & newTimes < t(2));


%% Plot random sample
old_rand = randperm(size(old_int,2),10);
old_plot = old_int(:,old_rand);

new_rand = randperm(size(new_int,2),10);
new_plot = new_int(:,new_rand);

%showSequences(ptOld,19,old_plot,[],0,[]);

%showSequences(ptNew,19,new_plot,[],0,[]);

%% compare new to random sample in some other time
t2 = [30*3600 31*3600];
new_int_other = newSeqs(:,newTimes > t2(1) & newTimes < t2(2));
new_plot_other = new_int_other(:,randperm(size(new_int_other,2),10));
%showSequences(ptNew,19,new_plot_other,[],0,[]);

old_int_other = oldSeqs(:,oldTimes > t2(1) & oldTimes < t2(2));
old_plot_other = old_int_other(:,randperm(size(old_int_other,2),10));
showSequences(ptOld,19,old_plot_other,[],0,[]);


end