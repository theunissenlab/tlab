function [amp_env] = make_amp_env(wrip_x, wrip_t, f_val, frameCount, amp_samp_rate, ind_low, ind_high, ranphase_flag);
% Generates the spectrogram for a single ripple given by wrip_x given in
% cycles/kHz and wrip_t given in Hz

nt = frameCount;                 % Number of points in time dimension
nbands = length(f_val);          % Number of points in frequency dimension

amp_env = zeros(nbands,nt); 

if (ranphase_flag) 
   phase_phi = 2*pi*rand; 
else 
   phase_phi =0; 
end 

t = 1:nt;
for tones=ind_low:ind_high   
    amp_env(tones,:) = cos(2*pi*wrip_x*f_val(tones)/1000.0 + 2*pi*wrip_t.*t/amp_samp_rate + phase_phi);    
end

%figure(3);
%imagesc(amp_env)
%axis xy;
%pause(3)
