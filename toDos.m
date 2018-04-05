%{

TO DOs:
- try longer times. Would be awesome to capture multiple seizures in one
plot!
- decide whether it makes sense to do a cleaning step (it's removing
dissimilar sequences, which is potentially removing an important part of
the data)


- [] validate the spike detectors
- [] remove EKG artifact?
- [] figure out why it runs out of PermGen memory on Borel
- [] confirm hard coded values - maybe ask Hoameng or Lohith what the voxel
to mm conversion is in the electrode location file?
- 

- by far the longest time period comes in grabbing data from ieeg.org. If I
can save gdf's that I have a decent amount of faith in, then this will be
very helpful.

Important info about the Tomlinson paper:
- about 4 spikes per second (about 100,000 spikes in 400 minutes)
- about 0.07 sequences per second (about 4 sequences per minute)
- so about every 60 spikes there is a sequence

- [] try source localization instead? I may be falsely assuming that
wherever I see the spike is where the actual source of the spike is.
Instead it could be that there is a spike below the surface that I am
seeing at multiple locations on the cortex at the same time

- [x] figure out why I'm mostly detecting EKG artifact-seems to be patient
specific??
     [x] make a script to more easily visualize detected spikes
     [x] try taking out the change i made to the Janca code to see if it is
     still a problem - yes
     [x] try it on a diff patient - not as bad on HUP80
     [x] write a script to take out spikes that occur too close to spikes
     detected on the ekg channel
- [x] try incorporating the Vanleer methods
- [x] correct some of the hard coded numbers to make sure physiologic. I
downloaded a paper Hirsh1991 (Synaptic physiology of horizontal connections
in th cat's visual cortex), which says that CV is about 0.3-1 m/s
- [x] should clean up the electrode localizer so that I am not remaking the
file each time; that's totally unnecessary (but then again it doesn't take
much time)
- [x] make it self contained so that it doesn't depend on outside tools (so
someone could download it and run it easily and so that I can run it on Borel)
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