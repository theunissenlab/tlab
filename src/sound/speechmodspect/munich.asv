% Tests filtering speech from the modulation spectrum

% Loads in speech
speech_name = 'C:\Documents and Settings\Frederic\My Documents\Matlab\speechmodspect\munich2';
[speech_in fs]= wavread(speech_name);

sound_out = modfilter(sound_in, fs, 32, 3, 0.0005, 500, 0, 0);
max_val = max(max(sound_out), abs(min(sound_out)));
sound_out = sound_out./max_val;
% Plays Filter Speech
soundsc(sound_out, fs);
speech_name_out = 'C:\Documents and Settings\Frederic\My Documents\Matlab\speechmodspect\munich_wf05';


wavwrite(sound_out, fs, 16, speech_name_out);
