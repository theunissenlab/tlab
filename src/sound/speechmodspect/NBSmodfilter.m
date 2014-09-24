% Tests filtering speech from the modulation spectrum

% Loads in speech
%song_name = 'C:\Documents and Settings\Frederic\My Documents\Data\birdsongs\ucsf\zfa_25\zfa_25.6.filt.wav';
song_name = 'C:\Users\Frederic\Documents\Data\Speech\Speech WavFiles\cheese_0.wav';

[song_in fs]= wavread(song_name);


soundsc(song_in, fs);

% Making SSI knock in
% High-pass sound above 20 Hz.

% First low pass filter below 5 Hz filter and 8 cyc/kHz
sound_lp = modfilter(song_in, fs, 32, 3, 0.008, 5, 0, 0);

% Then high-pass filter to above 0 Hz and 0.25 cyc/kHz
sound_hp = modfilter(sound_lp, fs, 32, 2, 0.00025, 0, 0, 0);

% Then prevent +1 and -1
sound_lp_fixed = sound_hp;
sound_lp_fixed(find(sound_lp<-1)) = -1;
sound_lp_fixed(find(sound_lp>1)) = 1;


soundsc(sound_lp_fixed, fs);

%song_name_out = 'C:\Documents and Settings\Frederic\My Documents\Data\birdsongs\ucsf\zfa_25\zfa_25_NBS_RAW.wav';
song_name_out = 'C:\Users\Frederic\Documents\Data\Speech\Filtered\Speech\cheese_NBSR.wav';

wavwrite(sound_lp_fixed, fs, 16, song_name_out);

% Then compress
m = max(sound_lp_fixed);
sound_lp_comp = fsig(sound_lp_fixed, m) - fsig(-sound_lp_fixed, m);

%song_name_out = 'C:\Documents and Settings\Frederic\My Documents\Data\birdsongs\ucsf\zfa_25\zfa_25_NBS_COMP.wav';

song_name_out = 'C:\Users\Frederic\Documents\Data\Speech\Filtered\Speech\cheese_NBSC.wav';

wavwrite(sound_lp_fixed, fs, 16, song_name_out);




