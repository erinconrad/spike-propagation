%{


To do:
- [] think about whether I should be using pre-ictal data immediately prior
to the seizure, more removed, or both
- [] make sure I calculated rms correctly
- [] I'm mainly detecting EKG artifact...
- [] make sure I did permutation test correctly

Vanleer et al. Millimeter-scale epileptiform spike propagation
patterns and their relationship to seizures. J Neural Eng. 2016; 13(2).

This paper looks at very small electrode arrays and looks at the
propagation patterns of spikes during seizures and outside of seizures.
They did this in cats.

What they did:
1) Detected spikes using their own method (I should probably not do this
since they are using a very different electrode array, animal, etc.).
2) On detection of a spike on any individual channel, take a 50 ms data
segment from all channels (2 ms prior to the amplitude crossing the spike
threshold and 48 ms post crossing). 
3) Calculate spike delay features for each segment: identify the relative
time of the maximum amplitude of the signal in the segment captured on each
channel with respect to the time of the maximum amplitude on the channel
with the earliest identifiable spike peak in the entire electrode group.
4) Also get power features by calculating root mean square values of each
channel after subtracting the window mean
5) Then independently normalize the delay maps and the power maps to force
it between 0 and 1.
6) Concatenate the 2 sets of features (360 delays, 360 powers), yielding
720 features per detected spike
7) Use PCA to retain the number of dimensions necessary to account for 99%
of the variance in the data
8) Do k-medians clustering (k=10) of the spike feature set to group spikes
(look at the paper for details)
9) to test that ictal spike patterns were different from interictal spike
patterns, they had a null hypothesis that the proportion of spikes occuring
during seizures was equal across all clusters. They used a permutation test
(n = 1 million), held cluster membership of each spike pattern fixed but
permuted the seizure/non-seizure labels for 1 million permutations,
recording the proportion of patterns with a seizure label within each
cluster
10) did a supervised analysis to predict interictal versus ictal with a QDA
classifier


%}