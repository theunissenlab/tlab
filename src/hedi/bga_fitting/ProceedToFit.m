function res_s=ProceedToFit(sp,freqs,freq_a,additional_guess)
% fit the spectrum using own guess and additional guess if available 
fmma=3500;
% 
guess=[];
guess=[additional_guess;guess];
pwr0=1;

[alpha,mu,sigma1,sigma2]=GuessFromSpectrum(sp,freqs,freq_a,fmma);
na=length(alpha);
for n=1:na;
    alpha_g=alpha(n);
    if ~isnan(alpha_g);               
        guess=[guess;alpha_g,mu,sigma1,sigma2,pwr0,1000];
    end    
end

[nguess,nparam]=size(guess);
res_s=[];

for i=1:nguess;
    gu=guess(i,:);    
    [alpha_f,mu_f,sigma_f1,sigma_f2,pwr_f,d]=GetMotoSpec(sp,gu(1),gu(2),gu(3),gu(4),gu(5));
    res_s=[res_s;alpha_f,mu_f,sigma_f1,sigma_f2,pwr_f,d];
end;
