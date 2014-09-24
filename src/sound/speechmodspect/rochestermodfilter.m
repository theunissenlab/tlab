% Tests filtering speech from the modulation spectrum

% Loads in speech
speech_name = 'C:\Documents and Settings\Frederic\My Documents\Data\Rochester\1minute.aif';
[speech_in fs]= aiffread(speech_name);
% Make space for output
sound_out = zeros(size(speech_in));

% Hardwired divisions for speech_name only
nsegs = 16;
segs = [ 0, ...
    119800, ...
    212400, ...
    315000, ...
    416500, ...
    520300, ...
    613000, ...
    719500, ...
    822000, ...
    918500, ...
    1021200, ...
    1113200, ...
    1211400, ...
    1320000, ...
    1412300, ...
    1478467];

for i=1:nsegs-1
    sound_in = speech_in(segs(i)+1:segs(i+1));
    % Plays Speech
    %soundsc(sound_in, fs);

    % Filter sound
    sound_out = modfilter(sound_in, fs, 128, 3, 0.0005, 500, 0, 0);

    % Plays Filter Speech
    %soundsc(sound_out, fs);

    % Stuff sound out
    speech_out(segs(i)+1:segs(i+1)) = sound_out;
    %pause;
end

speech_name_out = 'C:\Documents and Settings\Frederic\My Documents\Data\Rochester\1minute_wf05.aif';
aiffwrite(speech_out', fs, 16, speech_name_out);
