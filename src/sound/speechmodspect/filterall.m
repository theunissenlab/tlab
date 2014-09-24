% Perform a notch filter on all 100 sentences

soundfilename = '../WavFiles/writes_0.wav';
filename=soundfilename(1:(size(soundfilename,2)-6));
[sound_in, samprate, nbits]= wavread(soundfilename)


wf_high = 0.000;
wt_high = 0.000;
wf_it = 0.004;
wt_it = 4.000;

for tloop=1:5
    for floop=1:4

        wf_high=wf_high+wf_it;
        % Write sound files
        eval(['wavwrite(sound_out,' num2str(samprate) ',' '''' filename '_' filtertext '_t' num2str(tloop*4) 'f' num2str(floop*4) '.wav'')']);
    end

    wf_high=0.000;
    wt_high=wt_high+wt_it;

end