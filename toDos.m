%{

TO DOs:
- [] make it self contained so that it doesn't depend on outside tools (so
someone could download it and run it easily and so that I can run it on Borel)
- [] think about the value of these different cleaning measures, and think
about why I get such a different result when I restrict ties versus not (I
saved two sample versions in the results, one restricting ties and one not
restricting ties)
- [] decide about 2D vs 3D distance measure (right now everything is 3D)
- [] correct some of the hard coded numbers to make sure physiologic. I
downloaded a paper Hirsh1991 (Synaptic physiology of horizontal connections
in th cat's visual cortex), which says that CV is about 0.3-1 m/s
- [] try to validate the spike detectors: I need annotated data
- [] remove EKG artifact? Or remove detected spikes that are detected at
the same time in the EKG channel?
- [] should make it nicer for doing it for multiple patients
- [] should clean up the electrode localizer so that I am not remaking the
file each time; that's totally unnecessary (but then again it doesn't take
much time)
- [] try source localization instead? I may be falsely assuming that
wherever I see the spike is where the actual source of the spike is.
Instead it could be that there is a spike below the surface that I am
seeing at multiple locations on the cortex at the same time



- [x] nice way to plot recruitment latency maps
- [x] add in the potential to use some of Sam's stuff in my code, sequence
similarity matrix
- [x] add the ability to use either Sam or Hoameng's spike detector
- [x] add code to pull in seizure times and compare across hours.
- [x] come up with nice way to visualize sequences: this is a new script
called visualizeChLocations
- [x] find a way to visualize how similar the paths are: this is a script
called visualizeAvgPath
- [x] think about the fact that so many spikes seem to occur at the same
time. Think about the time steps and confirm that I am not losing sig figs.
It looks like the closest time steps are 0.005 seconds, so sampling at 200
Hz, which is what Sam's paper says.: The spike detector I was using
automatically downsamples the data to 200 Hz, even though we sample at 512.
I edited this code and then I also edited it so that we are using a
different filter algorithm to make the detector work better
- [x] can I optimize spike detection algorithm for the higher sampling
rate?: I just did this by switching the option for which filter to use


%}