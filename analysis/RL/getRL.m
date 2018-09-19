function [RL_avg,RL_each] = getRL(seqs)

RL_each = seqs - min(seqs,[],1);
RL_avg = nanmean(RL_each,2);

end