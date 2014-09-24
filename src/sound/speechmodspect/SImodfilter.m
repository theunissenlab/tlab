% Tests filtering speech from the modulation spectrum

% Loads in speech
song_name = 'C:\Documents and Settings\Frederic\My Documents\Data\birdsongs\ucsf\zfa_25\zfa_25.6.filt.wav';
[song_in fs]= wavread(song_name);


soundsc(song_in, fs);

% Making SSI knock in
% High-pass sound above 20 Hz.

% First low pass filter below 10 Hz filter and 0.25 cyc/kHz
sound_lp = modfilter(song_in, fs, 32, 3, 0.00025, 10, 0, 0);

% Then prevent +1 and -1
sound_lp_fixed = sound_lp;
sound_lp_fixed(find(sound_lp<-1)) = -1;
sound_lp_fixed(find(sound_lp>1)) = 1;
soundsc(sound_lp_fixed, fs);

song_name_out = 'C:\Documents and Settings\Frederic\My Documents\Data\birdsongs\ucsf\zfa_25\zfa_25_SI_RAW.wav';
wavwrite(sound_lp_fixed, fs, 16, song_name_out);
