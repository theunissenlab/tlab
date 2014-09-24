% Perform sound -> spectrogram -> modspect -> sound
% i.e. bypass filtering phase to see what distortion the 
% spectrogram inversion routine adds

function noFilter(sent)

disp(strcat('Performing spectrum unfiltered spectrum inversion on the ', strcat(' ', sent), ' sentence')); 
cd /auto/fdata/deep/SpeechTest-Syllables

wf_high = 0.000;
wt_high = 0.000;
wf_it = 0.004;
wt_it = 4.000;

soundfilename = strcat(sent, '/', sent, '_0', '.wav');
filename=soundfilename(1:(size(sent,2))+1);
[sound_in, samprate, nbits]= wavread(soundfilename);

filtertext = '_uf'; 

cd /auto/fhome/deep/speechmodspect
which modfilter
sound_out = modfilter(sound_in, samprate, 3, wf_high, wt_high, wf_it, wt_it);

cd /auto/fdata/deep/SpeechTest-Syllables
eval(['wavwrite(sound_out,' num2str(samprate) ',' '''' strcat(strrep(filename, '/WavFiles', '/SpeechTest-Syllables'), '/', sent) filtertext '.wav'')']);

wname = strrep(filename, '/', '');
wname = strcat(wname, filtertext, '.000');
fid = fopen(strcat('/auto/fdata/deep/SpeechTest-Syllables/', sent, '/', wname), 'w+');
binSoundOut = sound_out.*2^15;
binSoundOut = fix(binSoundOut);
fwrite(fid, binSoundOut, 'integer*2');
fclose(fid);





