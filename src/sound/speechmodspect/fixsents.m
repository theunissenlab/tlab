% Decrease power in all test sentences such that
% no clipping occurs after filtering

cd /auto/fhome/deep/SpeechTest-Sentences/
fid = fopen('Iowa.txt');
sents = textscan(fid, '%s%*[^\n]');
sents = sents{1};
sents = sort(sents);


for k = 1:length(sents)

    disp(strcat('Running ', num2str(k), 'of 100 trials')); 
    tic;
    soundfilename = strcat('../WavFiles/', sents{k}, '_0', '.wav');
    filename=soundfilename(1:(size(soundfilename,2)-6));
    [sound_in, samprate, nbits]= wavread(soundfilename);
    sound_in = sound_in./5; 
    wavwrite(sound_in, samprate, nbits, soundfilename); 
    
end
