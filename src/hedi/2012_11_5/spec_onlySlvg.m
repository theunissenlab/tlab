function [to, fo, logB] = spec_onlySlvg(sound_in, fband, samprate, dBScale)
% The Theunissen Lab Spectrogram with Guassian window.  Plots the
% spectrogram and oscillogram of the sound
% !!!! ATTENTION increment window fixée à 8 ms le 8juin2012 pour analyse DFA Spectros !!!!!!!!!!!!!!!!
% Ancienne valeur : 0.005 s

% Parameters for the Spectrogram
nstd = 6;
twindow = 1000*nstd/(fband*2.0*pi);           % Window length in ms - 6 times the standard dev of the gaussian window
winLength = fix(twindow*samprate/1000.0);  % Window length in number of points
winLength = fix(winLength/2)*2;            % Enforce even window length
increment = fix(0.008*samprate);           % Sampling rate of spectrogram in number of points - set at 500 Hz
f_low=0;                                 % Lower frequency bounds to get average amplitude in spectrogram
f_high=4000;                               % Upper frequency bound to get average amplitude in spectrogram
DBNOISE = dBScale;                          % dB in Noise for the log compression - values below will be set to zero.

% Calculate and plot the spectrogram    
[s, to, fo, pg] = GaussianSpectrum(sound_in, increment, winLength, samprate); 
logB = 20*log10(abs(s));
maxB = max(max(logB));
minB = maxB-DBNOISE;   

logB(find(logB<minB))=minB;     % replace every -Inf values by minB

imagesc(to*1000,fo,logB);          % to is in seconds
axis xy;
caxis('manual');
caxis([minB maxB]); 
cmap = spec_cmap();
colormap(cmap);

v_axis = axis; 
v_axis(3)=f_low; 
%v_axis(4)=f_high;
v_axis(4)=8000; % !!! Modifié le 10 mai 2012. Valeur par défaut = ligne du dessus.
axis(v_axis);                                

xlabel('Time (ms)'), ylabel('Frequency (Hz)');

end