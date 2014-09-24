function [to, fo, logB, pg, tError, fError] = spec(sound_in, fband, samprate, dBScale)
% The Theunissen Lab Spectrogram with Guassian window.  Plots the
% spectrogram and oscillogram of the sound

% Plot the oscillogram
soundlen = length(sound_in);
subplot(3,1,1);
t=1:soundlen;
t = (t-1).*1000/samprate;
plot(t,sound_in);
xlabel('time (ms)');
s_axis = axis;

% Plot the spectrogram
subplot(3,1,[2 3]);  
[to, fo, logB, pg, tError, fError] = spec_only(sound_in, fband, samprate, dBScale);

% Match the temporal axis for both plots
v_axis = axis; 
v_axis(1) = s_axis(1);
v_axis(2) = s_axis(2);
axis(v_axis);   

end