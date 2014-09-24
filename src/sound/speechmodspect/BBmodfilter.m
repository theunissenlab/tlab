% Tests filtering speech from the modulation spectrum

% Loads in speech
%song_name = 'C:\Documents and Settings\Frederic\My Documents\Data\birdsongs\ucsf\zfa_25\zfa_25.6.filt.wav';
song_name = 'C:\Users\Frederic\Documents\Data\Speech\Speech WavFiles\cheese_0.wav';

[song_in fs]= wavread(song_name);

% Make broadband noise with equal variance
%noise_in = normrnd(0, std(speech_in), length(speech_in), 1);

%sound_in = speech_in + noise_in;

soundsc(song_in, fs);

% Making FBOn knock in
% High-pass sound above 20 Hz.

% First high-pass filter to above 20 Hz and 0 cyc/kHz
sound_hp = modfilter(song_in, fs, 500, 2, 0, 20, 0, 0);

% Then prevent +1 and -1
sound_hp_fixed = sound_hp;
sound_hp_fixed(find(sound_hp<-1)) = -1;
sound_hp_fixed(find(sound_hp>1)) = 1;
soundsc(sound_hp_fixed, fs);

%song_name_out = 'C:\Documents and Settings\Frederic\My Documents\Data\birdsongs\ucsf\zfa_25\zfa_25_FBon_RAW.wav';

song_name_out = 'C:\Users\Frederic\Documents\Data\Speech\Filtered\Speech\cheese_FBBR.wav';

wavwrite(sound_hp_fixed, fs, 16, song_name_out);

% Then threshold
sound_hp_th = zeros(size(sound_hp_fixed));
ind_sound = find(sound_hp_fixed>0.5);

npts = length(ind_sound);
nint = fix(0.025*fs);
for i=1:npts
    if ( (sound_hp_fixed(ind_sound(i)-1) < sound_hp_fixed(ind_sound(i))) && ...
            (sound_hp_fixed(ind_sound(i)+1) < sound_hp_fixed(ind_sound(i))) )
        % We found a peak
        if ( ind_sound(i) - nint > 0)
            sound_hp_th(ind_sound(i)-nint:ind_sound(i)) = sound_hp_fixed(ind_sound(i)-nint:ind_sound(i));
        else
            sound_hp_th(1:ind_sound(i)) = sound_hp_fixed(1-nint:ind_sound(i)); 
        end
        if ( ind_sound(i) + nint < length(sound_hp_fixed) )
            sound_hp_th(ind_sound(i):ind_sound(i)+nint) = sound_hp_fixed(ind_sound(i):ind_sound(i)+nint);
        else
            sound_hp_th(ind_sound(i):end) = sound_hp_fixed(ind_sound(i):end); 
        end
    end
end
%song_name_out = 'C:\Documents and Settings\Frederic\My Documents\Data\birdsongs\ucsf\zfa_25\zfa_25_FBon_TH.wav';
song_name_out = 'C:\Users\Frederic\Documents\Data\Speech\Filtered\Speech\cheese_FBBT.wav'
wavwrite(sound_hp_th, fs, 16, song_name_out);

          
