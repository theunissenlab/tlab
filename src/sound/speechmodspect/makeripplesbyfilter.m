% Tests filtering speech from the modulation spectrum

% Generate gaussian white noise

tdur = 2;   % 2s sound
fs = 44100; % CD sampling rate
npts = fix(fs*tdur);
fhigh = 0.002;   % 2 cycles per kHz or 1 cyc/kHz
thigh = 30;       % 5 Hz or 30
nstim = 2;
tone_ramp = 500;   % 500 ms ramp
fnameRoot = '/Users/frederictheunissen/Documents/Data/Ripple Sounds/highhighRipple';

pd = makedist('Normal');

for i = 1:nstim
    
    soundNoise = random(pd, [1, npts]);
    soundNoise = addramp(soundNoise, fs, tone_ramp);

    soundsc(soundNoise, fs);
    
    pause(1.0);

    % Low pass filter the sound
    % soundNoiseFilt = songfilt(soundNoise, fs, 500, 10000, 0);
    % soundsc(soundNoiseFilt, fs)
    
    
    % First high pass filter to get rid of low spectral modulations
    %soundRippleTemp = modfilter(soundNoiseFilt, fs, 125, 2, flow, tlow, 0, 0);
    
    % soundsc(soundRippleTemp, fs);

    % First filter below 5 Hz filter and 0.1 cyc/kHz
    soundRipple = modfilter(soundNoise, fs, 125, 3, fhigh, thigh, 0, 0);

    % Filter and play back   
    soundRipple = songfilt(soundRipple, fs, 250, 10000, 0);
    soundsc(soundRipple, fs);
    audiowrite(sprintf('%s%d.wav', fnameRoot, i), soundRipple, fs);
    pause();
end

% song_name_out = 'C:\Documents and Settings\Frederic\My Documents\Data\birdsongs\ucsf\zfa_25\zfa_25_SI_RAW.wav';
% wavwrite(sound_lp_fixed, fs, 16, song_name_out);
