PittWavdir = '/auto/fdata/pgill/PittSounds/PittSounds';
%condir = '/auto/fdata/pgill/September_Cond/all_cond_stims/Spects';
%mlndir = '/auto/fdata/pgill/September_Cond/all_cond_stims/flatrip_Spects';
%%% Load Con Stims Here


%%% Load Pitt Stims Here

dirout = dir(fullfile(PittWavdir,'*.wav'));
short_list = [2    19    14     3    12     6     4    10     5    17    21     8     1    16     9 ...
    15    13    36    18    11     7    37    38    31    32    35    43    40    41    48 ...
    23    39    25    24    47    20];

spects_dir = '/auto/fdata/pgill/Zador/PittSpects';
fiatdir(spects_dir);
for jj = 1:length(short_list)
    [thisStim,fs] = wavread(fullfile(PittWavdir,dirout(short_list(jj)).name));
    stimResampled = resample(thisStim,97656,fs);
    stimResampled = stimResampled ./std(stimResampled);
    outstim = wav2spectrogramVarFs(stimResampled,97656,250);
    fname = fullfile(spects_dir,strrep(dirout(short_list(jj)).name,'.wav','.mat'));
    save(fname,'outstim');
end