August 27 2009

This readme file composed by Patrick Gill (patrick.robert.gill@gmail.com), currently at Cornell University.

To see an example surprise calculation, check out the code in "Demo_surprise_code.m".

In short, you need to create one MATLAB cell array of your spectrograms "stims" and a second cell array of stimulus names called "stim_names".

Then, you need to specify which stimuli you want outputs for with the array of 1s and 0s called "stim_used".  
If you want to know the surprise of every sound you entered, just say "stim_used = ones(size(stim_names));".

Lastly, if you want to specify how familiar each stim is (i.e. if you want to tweak the relative contribution of a stimulus to the P(S,D) model), you can do that with "stim_familiarity".
If all stimuli are equally familiar, just type "stim_familiarity = ones(size(stim_names));".
If there are stimuli that you don't want to contribute at all to P(S,D) (e.g. you might want to create a model that's completely unfamiliar with an artificial stimulus) set the corresponding value of stim_familiarity to 0.

Then, to calculate each frequency band's surprise, run (just as in "Demo_surprise_code.m"):

out = faster_surprise_of_one_frequency_band_7D(f_band,stims,stim_names,stim_familiarity,stim_used);

You might want to alter the default size of the domain used, in which case you would type:

out = faster_surprise_of_one_frequency_band_7D(f_band,stims,stim_names,stim_familiarity,stim_used,w,d,b);

Here, "w" is the width of the domain (in spectrogram time bins), default = 3.
"d" is the minimum latency between the stimulus whose surprise you're evaluating and the start of the domain, default = 3.
"b" is the half-height (minus 1/2) of the domain in spectrogram frequency bins; default = 4 (which gives a domain a full 9 frequency bins tall, centered on the frequency whose surprise you're evaluating.

The first band (for the lowest frequency) you can use is with f_band=1, but this will correspond to features at the frequency band (1+b) where "b" is the half-height of the domain you're using (measured in frequency bands).
The highest f_band you can use is N-(2*b), which gives features at the frequency of the (N-b)th spectrogram band.

There are 2*b fewer frequency bands in the surprise version than in the spectrograms because the surprise code does not examine frequency bins where the domain hangs below or peaks above the spectrogram.  
(If you really want every frequency bin, put "b" rows of 0s above and below your spectrograms; but be aware that the statistics of the bottom few rows in your surprise representation might have slightly reduced surprise.)

Please feel free to use this code for scientific purposes, but let me know about any other uses you might have before going ahead.  This code is a re-packaged version of the surprise code which went into Gill et al. 2008, but I've added the ability to make different stimuli have different weights in the corpus (see "stim_familiarity"). 