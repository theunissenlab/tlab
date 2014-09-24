% Perform the highpass filter on all 100 sentences

%Get array of all sentence names
cd /auto/fdata/deep/SpeechTest-Sentences/
fid = fopen('Iowa.txt');
sents = textscan(fid, '%s%*[^\n]');
sents = sents{1};
sents = sort(sents);

for k = 1:1

    soundfilename = strcat(sents{k} , '/', sents{k}, '_0', '.wav');
    filename=soundfilename(1:(size(soundfilename,2)-6));
    [sound_in, samprate, nbits]= wavread(soundfilename);
    

    wf_high = 0.000;
    wt_high = 0.000;
    wf_it = 0.004;
    wt_it = 4.000;
    
    filtertext='hpf';
    tloop = 1;
    wf_high=wf_high+wf_it;
     
    
    for floop=1:3
        % filter sound
         cd /auto/fhome/deep/speechmodspect
        sound_out = modfilter(sound_in, samprate, 2, wf_high, wt_high);
        num2str(floop*4)
        num2str(20)
        % Write sound_out in .wav format
%         eval(['wavwrite(sound_out,' num2str(samprate) ',' '''' strcat(strrep(filename, '/WavFiles', '/SpeechTest-Sentences'), '/', sents{k}) '_' filtertext '_t' num2str(floop*4) 'f' num2str(20) '.wav'')']);
%         
%        
%         %Write sound_out in .000 format
%         wname = strrep(filename, '../WavFiles/', '');
%         wname = strcat(wname, '_' ,filtertext, '_t', num2str(floop*4), 'f', num2str(20), '.000');
%         fid = fopen(strcat('/auto/fhome/deep/SpeechTest-Sentences/', sents{k}, '/', wname), 'w+');
%         binSoundOut = sound_out.*2^15;
%         binSoundOut = fix(binSoundOut);
%         fwrite(fid, binSoundOut, 'integer*2');
%         fclose(fid); 
%         
         wf_high=wf_high+wf_it;
        %pause;
    end

    wt_high=wt_high+wt_it;
    wf_high= 0;
    
    for tloop=1:4
        % filter sound
         cd /auto/fhome/deep/speechmodspect
        sound_out = modfilter(sound_in, samprate, 2, wf_high, wt_high);
        num2str(16)
        num2str(tloop*4)
%         %Write sound_out in .wav format
%         eval(['wavwrite(sound_out,' num2str(samprate) ',' '''' strcat(strrep(filename, '/WavFiles', '/SpeechTest-Sentences'), '/', sents{k}) '_' filtertext '_t' num2str(16) 'f' num2str(tloop*4) '.wav'')']);
% 
%         %Write sound_out in .000 format
%         wname = strrep(filename, '../WavFiles/', '');
%         wname = strcat(wname, '_' ,filtertext, '_t', num2str(16), 'f', num2str(tloop*4), '.000');
%         fid = fopen(strcat('/auto/fhome/deep/SpeechTest-Sentences/', sents{k}, '/', wname), 'w+');
%         binSoundOut = sound_out.*2^15;
%         binSoundOut = fix(binSoundOut);
%         fwrite(fid, binSoundOut, 'integer*2');
%         fclose(fid); 
        
        wt_high=wt_high+wt_it;
        %pause;
    end
    
end