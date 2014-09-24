function fit_calls(fname);

[paths,name,ext]=fileparts(fname);
addpath ../sound;
addpath ../2012_11_5;

alpha0=-0.4;
[freq_a,slope,intercept]=FreqAsAlpha(alpha0);

fprintf('doing %s\n',fname);
song=fitWave(fname,0.1,freq_a);

fprintf('done...\n');
outdir='/auto/k6/hsoula/stress/new_output/';

outname=sprintf('fit.%s.mat',name);
save([outdir,outname],'song');

