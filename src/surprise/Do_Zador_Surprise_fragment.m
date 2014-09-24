%  Demo surprise code, Feb 8 2008, by Patrick Gill.
%
%  Please feel free to use this code for scientific purposes, but let me
%  know about any other uses you might have before going ahead.  This code
%  is a re-packaged version of the surprise code which went into Gill et
%  al. 2008 (still in the publication tubes as of today), but I've added
%  the ability to make different stimuli have different weights in the
%  corpus (see "stim_familiarity").  I've run the sample code below and it
%  seems to work as expected, but it's always a good idea to do a reality
%  check before going too far.
%
%  This code loads three spectrograms of zebra finch
%  song, then puts them into the correct data format for running a surprise
%  calculation on frequency band number 22.  With three spectrograms of
%  about 2000 samples each (2 sec long), this takes about 1 minute on a
%  fairly new machine with more than enough ram to hold everything in
%  memory.  The resulting variable "out" will be a struct, with the
%  stimulus names of the first and third stimuli (which here were the only
%  ones used in presentations, so we don't need to find how surprising
%  the features of the second stim are).  The stimulus familiarity numbers are how familiar I
%  think each of the stimuli are: the middle one is 20 times as familiar as
%  the last one.  Values here should be nonnegative, but 0 is OK.
%
%  By default, this version computes a 7-dimensional P(S,D) and P(S|D), as
%  opposed to the 10-dimensional one used in Gill et al.  To run the full
%  10-dimensional version, use the function
%  "surprise_of_one_frequency_band" instead of
%  "surprise_of_one_frequency_band_7D".  (The extra three dimensions held a
%  tiny bit in the Gill et al. CLM dataset, but they increase the time it
%  takes to compute surprise by a large factor.)
%
%  To do a whole spectrogram, you'll have to loop through each f_band
%  between the lowest (1) and the highest (end-(b*2)) where b is the number
%  of frequency bins in the domain (default = 4).

%{
loaded = load('008B5A5C0C3E76BFBF97278658FB6309.mat');
stims = {loaded.outSpectrum};
stim_names = {'008B5A5C0C3E76BFBF97278658FB6309.mat'};

loaded = load('1470489635dd93410408ce9f8fb2f7d9.mat');
stims{end+1} = loaded.outSpectrum;
stim_names{end+1} = '1470489635dd93410408ce9f8fb2f7d9.mat';

loaded = load('A10F5407D5D86F179E9305D1453632F7.mat');
stims{end+1} = loaded.outSpectrum;
stim_names{end+1} = 'A10F5407D5D86F179E9305D1453632F7.mat';

stim_familiarity = [4 10 .5];

stim_used = [1 0 1];
%}
stims = {};
stim_names = {};
rootSpecsDir = '/auto/fdata/pgill/Zador/Spectrograms_cont/250';
savedir = '/auto/fdata/pgill/Zador/Surprise/in_progress_cont';
fiatdir(savedir);
%dirout = dir(fullfile(rootSpecsDir,'class*'));

%for jj = 1:length(dirout)
    allSpecs = dir(fullfile(rootSpecsDir,'*.mat'));
    for kk =1:length(allSpecs)
        loaded = load(fullfile(rootSpecsDir,allSpecs(kk).name));
        stims{end+1} = loaded.outstim;
        stim_names{end+1} = [allSpecs(kk).name];
    end
%end
stim_familiarity = ones(1,length(stims));
stim_used = stim_familiarity;
for f_band = 1:86%94
    %f_band = 22;
    try
        savename = fullfile(savedir,[num2str(f_band) '.mat']);
        if ~exist(savename,'file')
            out = surprise_of_one_frequency_band_7D(f_band,stims,stim_names,stim_familiarity,stim_used);
            save(savename,'out');
        end
    catch
        disp(lasterr);
    end
end
%  So, if you wanted 10-D characterizations of P(S,D) and P(S|D), you'd use
%  the following instead:
%  out = surprise_of_one_frequency_band(f_band,stims,stim_names,stim_familiarity,stim_used);
