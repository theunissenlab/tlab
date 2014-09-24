function [info, noiseentropy, totalentropy]=make_psth_gammanoonset(spiketrain);
%[info, noiseentropy, totalentropy]=make_psth_gammanoonset(spiketrain);
%returns 5 values of info, noiseentropy, totalentropy, at extrap, 16, 12, 8 4 ms wordlengths.
%this smooths each spike train with a varying-window-length gaussian.
[meanrate] = gauss_filter_vwbytrial(spiketrain,  5);
meanrate=mean(meanrate);
%meanrate = gauss_filter_psth_near(mean(spiketrain),  1);
%the following outputs the estimated meanrate into a file called rate, to
%be called by noisyspike.

%talfitdata is the time rescaled ISIs of spiketrain

  [talvecnew, talvecold, rescaled_spikes]=time_rescaled_tals_JN(spiketrain, 5);
%gamfit gets the estimated constant gamma that best describes the ISI
%distributions
if or(isempty(talvecnew), length(talvecnew)<5)
    %assume poisson if there's not enough data
    gammaconst=1;
    gammaerror=1;
    info=zeros(5,1);
noiseentropy=zeros(5,1);
totalentropy=zeros(5,1);
else

[paramdata pci]=gamfit(talvecnew);
%gammaerrror is the percentage of estimate that is at the 5% uncertaintly
%level 

gammaconst=paramdata(1);
spikyneuron1=make_spikes3(meanrate, gammaconst, 500);
info=zeros(4,1);
noiseentropy=zeros(4,1);
totalentropy=zeros(4,1);
for k = 4:-1:1
    wordlength=2*k;
    [infot,totalentropyt, noiseentropyt]= get_info(spikyneuron1, wordlength, 2, 2);
    info(k)=infot;
    noiseentropy(k)=noiseentropyt;
    totalentropy(k)=totalentropyt;
end
ifit=polyfit(1./[16 12 8 4], info(1:4)', 1);
tfit=polyfit(1./[16 12 8 4], totalentropy(1:4)', 1);
nfit=polyfit(1./[16 12 8 4], noiseentropy(1:4)', 1);

totalentropy=[polyval(tfit, 0); totalentropy];
noiseentropy=[polyval(nfit, 0); noiseentropy];
info=[totalentropy(1)-noiseentropy(1); info];
end
