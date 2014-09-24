function sound_out = addramp(sound_in, samp_rate, ramp_duration)
% Adds a cosine onset and offset ramp to the sound of length ramp_duration
% in ms

samp_ms = samp_rate/1000;
sound_out = sound_in;
lensound = length(sound_in);

nend = fix(ramp_duration*samp_ms);
if (nend >= lensound)
    nend = lensound -1;
end

% Onset ramp
for t = 1:nend; 
   mult1=0.5*(1.0-cos(pi*t/nend)); 
   sound_out(t) = sound_out(t)*mult1; 
end 

nbeg = lensound - nend;
for t = nbeg : lensound 
   mult1=0.5*(1.0-cos(pi*(lensound-t)/nend)); 
   sound_out(t) = sound_out(t)*mult1; 
end 