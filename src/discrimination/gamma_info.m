function [info, noiseentropy, totalentropy, gammaconst, bandwidth, rategamma, avgff]=gamma_info(spiketrain)
% Calculates the gamma information as well as the bandwidth and the entropy
% of the mean rate obtained by convolving the psth with a varying guassian
% window.
% 5 values of info, noiseentropy, totalentropy are returned for: extrap, 16, 12, 8 4 ms wordlengths.

%this smooths each spike train with a varying-window-length gaussian.
[meanrate] = gauss_filter_vwbytrial(spiketrain,  5);
meanrate=mean(meanrate);

% The bandwidth of the meanrate function
bandwidth=make_bandwidth(meanrate);

% The Fano factor calculated with 30 ms window
[me, sd, v] = cv2(spiketrain, 30);
ff=v(find(me))./me(find(me));
avgff=mean(ff);

% The gamma factor for the rate distribution
[gamtemp ratepci]=gamfit(meanrate);
%rategammaconst is the fit of the pdf of meanrate function.
%the bigger the rategammaconst, the more normal
rategamma=gamtemp(1);
% rategammalb=ratepci(1,1);
% rategammaub=ratepci(2,1);
% rategammaerror=(rategammaub-rategamma)/rategamma;


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
% gammaconstlb=pci(1,1);
% gammaconstub=pci(2,1);
% gammaerror=(gammaconstub-gammaconst)/gammaconst;

% Generate Random data
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
%Modified by Jonathan Shih 2009-10-12 . . . extrapolation had incorrect
%order of word lengths
% ifit=polyfit(1./[16 12 8 4], info(1:4)', 1);
% tfit=polyfit(1./[16 12 8 4], totalentropy(1:4)', 1);
% nfit=polyfit(1./[16 12 8 4], noiseentropy(1:4)', 1);
ifit=polyfit(1./[4 8 12 16], info(1:4)', 1);
tfit=polyfit(1./[4 8 12 16], totalentropy(1:4)', 1);
nfit=polyfit(1./[4 8 12 16], noiseentropy(1:4)', 1);

totalentropy=[polyval(tfit, 0); totalentropy];
noiseentropy=[polyval(nfit, 0); noiseentropy];
info=[totalentropy(1)-noiseentropy(1); info];


end
