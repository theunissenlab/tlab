cd /auto/fhome/deep/SpeechTest-Sentences/
fid = fopen('Iowa.txt');
sents = textscan(fid, '%s%*[^\n]');
sents = sents{1};
sents = sort(sents);

for k = 1:44
    cdcom = strcat('cd /auto/fhome/deep/SpeechTest-Sentences/', sents{k}); 
    eval(cdcom); 
    
    name1 = strcat(sents{k}, '_hpf_t8f12.wav');
    delete(name1); 
    
    name2 = strcat(sents{k}, '_hpf_t4f12.wav');
    delete(name2); 
    
    name3 = strcat(sents{k}, '_hpf_t4f8.wav');
    delete(name3); 
    
    name4 = strcat(sents{k}, '_hpf_t4f4.wav');
    delete(name4); 
    
    name5 = strcat(sents{k}, '_hpf_t12f12.wav');
    delete(name5); 
    
    name6 = strcat(sents{k}, '_nf_t20f16.wav');
    delete(name6); 
    
    name7 = strcat(sents{k}, '_nf_t20f12.wav');
    delete(name7); 
    
    name8 = strcat(sents{k}, '_nf_t20f8.wav');
    delete(name8); 
    
    name9 = strcat(sents{k}, '_nf_t20f4.wav');
    delete(name9); 
    
end
