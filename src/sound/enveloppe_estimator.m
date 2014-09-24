function amp_env = enveloppe_estimator(sound_in, fs, f_high, fs_env)
% Calculates the amplitude envelope of sound_in, sampled at fs by
% rectifying and low-pass filtering with a cutt-off frequency of fhigh. The
% amplitude enveloppe is then resampled at famps

% Rectify
sound_rect = abs(sound_in);

% Low pass filter
% Figure out length of filter
nframes = length(sound_in);
if ( nframes > 3*512 ) 
    nfilt = 512;
elseif ( nframes > 3*64 ) 
    nfilt = 64;
    
elseif ( nframes > 3*16 ) 
    nfilt = 16;
else
    error('Song data is too short for filtering');
end

% Generate filter and filter rectified enveloppe
lowpass_filter = fir1(nfilt, f_high*2.0/fs);
amp_env_full = filtfilt(lowpass_filter, 1, sound_rect);

% Resample to desired sampling rate
if fs ~= fs_env
    amp_env = resample(amp_env_full, fs_env, fs);
else
    amp_env = amp_env_full;
end

return