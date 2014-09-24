%% Read all the wav files in the directory and make spectrograms

wavfiles = dir('*.wav');

nfiles = length(wavfiles);

ncols = 4;
nrows = ceil(nfiles/ncols);
fband = 100;
dBScale = 50;
figure();

Ptot = zeros(129, 1);

for i=1:nfiles
    fileName = wavfiles(i).name;
    [sound fs] = wavread(fileName);
    
    subplot(nrows, ncols, i);
    spec_only(sound, fband, fs, dBScale);
    [Pxx, f] = pwelch(sound,128,64,256,fs,'onesided');
    Ptot = Ptot + Pxx;
end

Ptot = Ptot./nfiles;