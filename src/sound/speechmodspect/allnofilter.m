% Perform sound -> spectrogram -> modspect -> sound
% i.e. bypass filtering phase to see what distortion the 
% spectrogram inversion routine adds

function noFilter(k)

cd /auto/fdata/deep/SpeechTest-Sentences/
fid = fopen('Iowa.txt');
sents = textscan(fid, '%s%*[^\n]');
sents = sents{1};
sents = sort(sents);


filtertext='_uf';

wf_high = 0.000;
wt_high = 0.000;
wf_it = 0.004;
wt_it = 4.000;

for k = 1:1

    soundfilename = strcat(sents{k} , '/', sents{k}, '_0', '.wav');
    filename=soundfilename(1:(size(soundfilename,2)-6));
    [sound_in, samprate, nbits]= wavread(soundfilename);
    
    cd /auto/fhome/deep/speechmodspect
    sound_out = modfilter(sound_in, samprate, 3, wf_high, wt_high, wf_it, wt_it);
    
%     eval(['wavwrite(sound_out,' num2str(samprate) ',' '''' strcat(strrep(filename, '/WavFiles', '/SpeechTest-Sentences'), '/', sents{k}) filtertext '.wav'')']);
% 
%     wname = strrep(filename, '../WavFiles/', '');
%     wname = strcat(wname, filtertext, '.000');
%     fid = fopen(strcat('/auto/fhome/deep/SpeechTest-Sentences/', sents{k}, '/', wname), 'w+');
%     binSoundOut = sound_out.*2^15;
%     binSoundOut = fix(binSoundOut);
%     fwrite(fid, binSoundOut, 'integer*2');
%     fclose(fid); 
    
end

    
    
    