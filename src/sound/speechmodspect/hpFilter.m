% Perform the highpass filter on all 100 sentences
function hpFilter(sent)

cd /auto/fdata/deep/SpeechTest-Syllables/
disp(strrep(strcat('Performing high pass filter on theX', strcat(' ', sent), ' syllable'), 'X', ' ')); 

soundfilename = strcat(sent, '/', sent, '_0', '.wav');
filename=soundfilename(1:(size(sent,2))+1);
[sound_in, samprate, nbits]= wavread(soundfilename);


wf_high = 0.000;
wt_high = 0.000;
wf_it = 0.004;
wt_it = 4.000;

filtertext='hpf';
tloop = 1;
wf_high=wf_high+wf_it;

num = 0; 

% for floop=1:3
%     % filter sound
%     cd /auto/fhome/deep/speechmodspect
%     which modfilter
%     sound_out = modfilter(sound_in, samprate, 2, wf_high, wt_high);
% 
%     cd /auto/fdata/deep/SpeechTest-Syllables
%     % Write sound_out in .wav format
%     eval(['wavwrite(sound_out,' num2str(samprate) ',' '''' strcat(strrep(filename, '/WavFiles', '/SpeechTest-Syllables'), '/', sent) '_' filtertext '_t' num2str(floop*4) 'f' num2str(20) '.wav'')']);
% 
% 
%     %Write sound_out in .000 format
%     wname = strrep(filename, '/', '');
%     wname = strcat(wname, '_' ,filtertext, '_t', num2str(floop*4), 'f', num2str(20), '.000');
%     fid = fopen(strcat('/auto/fdata/deep/SpeechTest-Syllables/', sent, '/', wname), 'w+');
%     binSoundOut = sound_out.*2^15;
%     binSoundOut = fix(binSoundOut);
%     fwrite(fid, binSoundOut, 'integer*2');
%     fclose(fid);
%       
%     wf_high=wf_high+wf_it;
%     
%     num = num + 1; 
%     
%     disp(strrep(strcat('FinishedX', strcat(' ', num2str(num)), ' of 7'), 'X', ' ')); 
%     
%     %pause;
% end

wt_high=wt_high+wt_it;
wf_high= 0;

for tloop=1:4
    % filter sound
    cd /auto/fhome/deep/speechmodspect
    which modfilter
    sound_out = modfilter(sound_in, samprate, 2, wf_high, wt_high);

    %Write sound_out in .wav format
    cd /auto/fdata/deep/SpeechTest-Syllables
    eval(['wavwrite(sound_out,' num2str(samprate) ',' '''' strcat(strrep(filename, '/WavFiles', '/SpeechTest-Syllables'), '/', sent) '_' filtertext '_t' num2str(16) 'f' num2str(tloop*4) '.wav'')']);

    %Write sound_out in .000 format
    wname = strrep(filename, '/', '');
    wname = strcat(wname, '_' ,filtertext, '_t', num2str(16), 'f', num2str(tloop*4), '.000');
    fid = fopen(strcat('/auto/fdata/deep/SpeechTest-Syllables/', sent, '/', wname), 'w+');
    binSoundOut = sound_out.*2^15;
    binSoundOut = fix(binSoundOut);
    fwrite(fid, binSoundOut, 'integer*2');
    fclose(fid);

    wt_high=wt_high+wt_it;
    
    num = num + 1; 
    disp(strrep(strcat('FinishedX', strcat(' ', num2str(num)), ' of 4'), 'X', ' ')); 
    %pause;
end

