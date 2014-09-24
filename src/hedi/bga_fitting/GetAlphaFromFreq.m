function alpha=GetAlphaFromFreq(freq_a,freq);

df=50; %Hz
d=freq_a(abs(freq_a(:,2)-freq)<df,1);
alpha=mean(d);

