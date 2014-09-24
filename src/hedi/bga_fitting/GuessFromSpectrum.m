function [alpha,mu,sigma1,sigma2]=GuessFromSpectrum(spectrum,sfreqs,freq_a,fmma)

[m,im]=max(spectrum);
fmax=sfreqs(im);
mu=fmax;
zmu=find(spectrum>m*0.15);
%sfreqs(zmu(end))
%sfreqs(zmu(1))
%sfreqs(zmu(end))-fmax;
%fmax-sfreqs(zmu(1));

sigma1=sfreqs(zmu(end))-fmax;sigma2=fmax-sfreqs(zmu(1));
alpha_m=GetAlphaFromFreq(freq_a,fmax);

% likely fundamental in the range [500;1000];
% look for ns inside this range 
nmax=fmax/500;


alpha=[alpha_m];
for n=2:nmax;
    if fmax/n<fmma;
        alpha=[alpha;GetAlphaFromFreq(freq_a,fmax/n)];
    end;
end;
