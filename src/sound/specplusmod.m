function [to, fo, logB] = specplusmod(sound_in, fband, samprate, dBScale)
% The Theunissen Lab Spectrogram with Guassian window.  Plots the
% spectrogram and oscillogram of the sound

% Parameters for the Spectrogram
nstd = 6;
twindow = 1000*nstd/(fband*2.0*pi);           % Window length in ms - 6 times the standard dev of the gaussian window
winLength = fix(twindow*samprate/1000.0);  % Window length in number of points
winLength = fix(winLength/2)*2;            % Enforce even window length
increment = fix(0.001*samprate);           % Sampling rate of spectrogram in number of points - set at 1 kHz
f_low=0;                                 % Lower frequency bounds to get average amplitude in spectrogram
f_high=5000;                               % Upper frequency bound to get average amplitude in spectrogram
DBNOISE = dBScale;                          % dB in Noise for the log compression - values below will be set to zero.


% Plot the oscillogram
figure();
soundlen = length(sound_in);
subplot(4,1,1);
t=1:soundlen;
t = (t-1).*1000/samprate;
plot(t,sound_in);
xlabel('time (ms)');
s_axis = axis;
s_axis(1) = 0;
s_axis(2) = t(end);
axis(s_axis);

% Calculate and plot the spectrogram
subplot(4,1,[2 3 4]);
[s, to, fo, pg] = GaussianSpectrum(sound_in, increment, winLength, samprate);
logB = 20*log10(abs(s));
maxB = max(max(logB));
minB = maxB-DBNOISE;
%pcolor(t.*1000,f,logB); shading interp;
imagesc(to*1000,fo,logB);          % to is in seconds
axis xy;
caxis('manual');
caxis([minB maxB]);
cmap = spec_cmap();
colormap(cmap);

%  Match the axis to oscillogram
v_axis = axis;
%v_axis(2)=4.0;
v_axis(1) = s_axis(1);
v_axis(2) = s_axis(2);
v_axis(3)=f_low;
v_axis(4)=f_high;
axis(v_axis);
xlabel('time (ms)'), ylabel('Frequency');
%axis off;

figure();
% Find the mean level of the amplitude envelope as a function of frequency
ampsrate = samprate/increment;
sabs = logB;
sabs(find(sabs <= minB)) = minB;
%Subtract DC
sabs = sabs - mean(mean(sabs));
nb = size(sabs,1);
nt = length(to);

% calculate the 2D fft
fabs = fft2(sabs);

% calculate amplitude and phase
amp_fabs = abs(fabs);
phase_fabs = angle(fabs);
amp_fabs_db = 20*log10(amp_fabs);
max_amp = max(max(amp_fabs_db));
min_amp = min(min(amp_fabs_db));

% Calculate the labels for temporal and spectral frequencies in physical
% units
% f_step is the separation between frequency bands
fstep = fo(2);
for ib=1:ceil((nb+1)/2)
    dwf(ib)= (ib-1)*(1/(fstep*nb));
    if (ib > 1)
        dwf(nb-ib+2)=-dwf(ib);
    end
end

% ampsrate is the sampling rate of the amplitude enveloppe
for it=1:ceil((nt+1)/2)
    dwt(it) = (it-1)*(ampsrate/nt);
    if (it > 1 )
        dwt(nt-it+2) = -dwt(it);
    end
end



imagesc(fftshift(dwt), fftshift(dwf).*1000, fftshift(amp_fabs_db));
caxis('manual');
caxis([max_amp-DBNOISE max_amp]);
title('Amplitude Spectrum');
axis xy;
axis([-60 60 0 20]);
colormap(cmap);

figure();
%imagesc(dwt, dwf.*1000, unwrap(unwrap(fftshift(phase_fabs)),1));
% Set phase to -3.15 for small amplitude
ind = find(amp_fabs_db < max_amp-DBNOISE);
phase_fabs(ind) = -3.15;

imagesc(fftshift(dwt), fftshift(dwf).*1000, fftshift(phase_fabs));
cmap = colormap('HSV');
cmap(1,:) = [1 1 1];
colormap(cmap);
colorbar();
title('Phase Spectrum');
axis xy;
axis([-60 60 0 20]);


 

end