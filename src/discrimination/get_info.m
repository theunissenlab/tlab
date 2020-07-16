function [info, totalentropy, noiseentropy]=get_info(spiketrain, wordlength, lettersize, rlettersize);
%[info, totalentropy, noiseentropy]=get_info(spiketrain, wordlength, lettersize, rlettersize);
%rlettersize is the number of bins that the spiketrain should be grouped together
%before continuing calculation.
ntrials=size(spiketrain,1);
ncolumns=size(spiketrain, 2); 
maxcolumns=ncolumns-mod(ncolumns,rlettersize);
spiketrain=spiketrain(:,1:maxcolumns);
ncolumns=size(spiketrain, 2);   
if rlettersize>1
    spiketrain2=reshape(spiketrain,ntrials, rlettersize,ncolumns/rlettersize);
    spiketrain2=sum(spiketrain2, 2);
    spiketrain2=reshape(spiketrain2, ntrials, ncolumns/rlettersize);
    spiketrain=spiketrain2;
    clear spiketrain2
end


 ncolumns=size(spiketrain, 2);   
maxcolumns=ncolumns-mod(ncolumns,wordlength);
spiketrain=spiketrain(:,1:maxcolumns);
totalentropy=entropy(spiketrain, wordlength);
noiseentropy=noise_entropy(spiketrain, wordlength, max(max(spiketrain))+1);

totalentropy=totalentropy/(wordlength*lettersize);
noiseentropy=noiseentropy/(wordlength*lettersize);
info=totalentropy-noiseentropy;
