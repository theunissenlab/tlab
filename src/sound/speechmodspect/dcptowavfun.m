function dcptowavfun(filename,ofilename)
[samplefreq nframes data_wave] = read_song(filename);

data_wave = data_wave/pow2(15); 
%   This effectively scales the dynamic range between -1 and 1 
%   and preserves the filtering we did in songfilt

wavwrite(data_wave,samplefreq,ofilename);