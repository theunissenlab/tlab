% Perform a notch filter on all 100 syllables
function notchFilter(sent)

cd /auto/fdata/deep/SpeechTest-Syllables/
disp(strrep(strcat('Performing notch filter on theX', strcat(' ', sent), ' sentence'), 'X', ' ')); 

soundfilename = strcat(sent, '/', sent, '_0', '.wav');
filename=soundfilename(1:(size(sent,2))+1);
[sound_in, samprate, nbits]= wavread(soundfilename);

cd /auto/fhome/deep/speechmodspect

filtertext='nf';

wf_high = 0.000;
wt_high = 0.000;
wf_it = 0.004;
wt_it = 4.000;

num = 0; 

for tloop=1:5
    for floop=1:4
        
        
        which modfilter
        sound_out = modfilter(sound_in, samprate, 1, wf_high, wt_high, wf_it, wt_it);
        
        wf_high=wf_high+wf_it;
        floop * 4
        tloop * 4
        
        cd /auto/fdata/deep/SpeechTest-Syllables
        
        % Write sound_out in .wav format
        eval(['wavwrite(sound_out,' num2str(samprate) ',' '''' strcat(strrep(filename, '/WavFiles', '/SpeechTest-Syllables'), '/', sent) '_' filtertext '_t' num2str(floop*4) 'f' num2str(tloop*4) '.wav'')']);

        %Write sound_out in binary format (16 bit signed integer)
        wname = strrep(filename, '/', '');
        wname = strcat(wname, '_' ,filtertext, '_t', num2str(floop*4), 'f', num2str(tloop*4), '.000');
        fid = fopen(strcat('/auto/fdata/deep/SpeechTest-Syllables/', sent, '/', wname), 'w+');
        binSoundOut = sound_out.*2^15;
        binSoundOut = fix(binSoundOut);
        fwrite(fid, binSoundOut, 'integer*2');
        fclose(fid);
        num = num + 1; 
        disp(strrep(strcat('FinishedX', strcat(' ', num2str(num)), ' of 20'), 'X', ' ')); 
     end

    wf_high=0.000;
    wt_high=wt_high+wt_it;

end
    


