cd /auto/fhome/deep/SpeechTest-Syllables
fid = fopen('syll.txt');
sents = textscan(fid, '%s'); 
sents = sents{1};
sents = sort(sents);

for k = 14:14
    noFilter(sents{k}); 
    hpFilter(sents{k}); 
    notchFilter(sents{k}); 
    
end