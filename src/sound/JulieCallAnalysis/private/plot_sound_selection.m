function plot_sound_selection(VocType,sound_in, samprate, on)

% Parameters for the Spectrogram
nstd = 6;
fband = 50;
twindow = 1000*nstd/(fband*2.0*pi);           % Window length in ms - 6 times the standard dev of the gaussian window
winLength = fix(twindow*samprate/1000.0);  % Window length in number of points
winLength = fix(winLength/2)*2;            % Enforce even window length
increment = fix(0.001*samprate);           % Sampling rate of spectrogram in number of points - set at 1 kHz
f_low=0;                                 % Lower frequency bounds to get average amplitude in spectrogram
f_high=8000;                               % Upper frequency bound to get average amplitude in spectrogram
DBNOISE = 40;                          % dB in Noise for the log compression - values below will be set to zero.



% Calculate and plot the spectrogram   
figure(1)
subplot(2,1,1); 
cla;
[s, to, fo, pg] = GaussianSpectrum(sound_in, increment, winLength, samprate); 
logB = 20*log10(abs(s));
maxB = max(max(logB));
minB = maxB-DBNOISE;            

imagesc(to,fo,logB);          % to is in seconds
axis xy;
caxis('manual');
caxis([minB maxB]); 

axis([0 length(sound_in)/samprate f_low f_high]);                                
ylabel('Frequency kHz');
cmap = spec_cmap();
colormap(cmap);


% Add a title  
subplot(2,1,1);
titleName = sprintf('Vocalization: %s', VocType);
title(titleName);

%plot waveline and on
subplot(2,1,2);
Max=0.1*ceil(10*max(sound_in));
plot(sound_in)
hold on
plot(Max*on, 'r');
hold on 
axis([0 length(on) -Max Max])
hold off;

