function [train,test] = makeTrainAndTest(ntotal,ntrain)

allIdx = [1:ntotal]';

train = randperm(ntotal,ntrain)';
test = find(~ismember(allIdx,train));

end