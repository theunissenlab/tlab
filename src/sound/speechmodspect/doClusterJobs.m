% Runs notchFilter, hpFilter, noFilter on all 100 sentences
% using 100 iterations for the spectrum inversion to write back the wav 
% files / binary files

cd /auto/fhome/deep/SpeechTest-Syllables
fid = fopen('syll.txt');
sents = textscan(fid, '%s'); 
sents = sents{1};
sents = sort(sents);

for k = 2:length(sents)
   
   notchcom = ['addpaths;notchFilter(''' sents{k} ''')']; 
   hpcom = ['addpaths;hpFilter(''' sents{k} ''')']; 
   nocom = ['addpaths;noFilter(''' sents{k} ''')'];
   
   dbaddqueuemaster(notchcom); 
   dbaddqueuemaster(hpcom);
   dbaddqueuemaster(nocom); 
   
end


