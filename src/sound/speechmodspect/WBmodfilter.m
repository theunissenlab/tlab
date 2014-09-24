% Tests filtering speech from the modulation spectrum

% Loads in speech
%song_name = 'C:\Documents and Settings\Frederic\My Documents\Data\birdsongs\ucsf\zfa_25\zfa_25.6.filt.wav';

song_name = 'C:\Users\Frederic\Documents\Data\Speech\Speech WavFiles\cheese_0.wav';
[song_in fs]= wavread(song_name);


soundsc(song_in, fs);

% Making SSI knock in
% First low pass filter below 30 Hz filter and 0.5 cyc/kHz
sound_lp = modfilter(song_in, fs, 250, 3, 0.0005, 30, 0, 0);

% Then high-pass filter to above 1 Hz and 0.1 cyc/kHz
% sound_hp = modfilter(sound_lp, fs, 250, 2, 0.0001, 5, 0, 0);

% Then prevent +1 and -1
sound_lp_fixed = sound_lp;
sound_lp_fixed(find(sound_lp<-1)) = -1;
sound_lp_fixed(find(sound_lp>1)) = 1;


soundsc(sound_lp_fixed, fs);

%song_name_out = 'C:\Documents and Settings\Frederic\My Documents\Data\birdsongs\ucsf\zfa_25\zfa_25_WB_RAW.wav';
song_name_out = 'C:\Users\Frederic\Documents\Data\Speech\Filtered\Speech\cheese_WBR.wav';


wavwrite(sound_lp_fixed, fs, 16, song_name_out);

% Then compress
m = max(sound_lp_fixed);
sound_lp_comp = fsig(sound_lp_fixed, m) - fsig(-sound_lp_fixed, m);

%song_name_out = 'C:\Documents and Settings\Frederic\My Documents\Data\birdsongs\ucsf\zfa_25\zfa_25_WB_COMP.wav';
song_name_out = 'C:\Users\Frederic\Documents\Data\Speech\Filtered\Speech\cheese_WBC.wav';


wavwrite(sound_lp_comp, fs, 16, song_name_out);
