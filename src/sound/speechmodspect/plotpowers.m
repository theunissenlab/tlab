function plotpowers

cd /auto/fdata/deep/SpeechTest-Sentences/
fid = fopen('Iowa.txt');
sents = textscan(fid, '%s%*[^\n]');
sents = sents{1};
sents = sort(sents);

cd /auto/fdata/deep/SpeechTest-Syllables/
fid = fopen('syll.txt');
sylls = textscan(fid, '%s%*[^\n]');
sylls = sylls{1};
sylls = sort(sylls);

powersents = []; 
powersylls = []; 

cd /auto/fdata/deep/SpeechTest-Sentences/
for k = 1:length(sents)

    soundfilename = strcat(sents{k} , '/', sents{k}, '_0', '.wav');
    filename=soundfilename(1:(size(soundfilename,2)-6));
    [sound_in, samprate, nbits]= wavread(soundfilename); 
    powersents = [powersents std(sound_in)]; 
end


cd /auto/fdata/deep/SpeechTest-Syllables/
for k = 1:length(sylls)

    soundfilename = strcat(sylls{k} , '/', sylls{k}, '_0', '.wav');
    filename=soundfilename(1:(size(soundfilename,2)-6));
    [sound_in, samprate, nbits]= wavread(soundfilename); 
    powersylls = [powersylls std(sound_in)]; 
end

plot(powersents);title('Sentences')
figure;plot(powersylls);title('Syllables')

