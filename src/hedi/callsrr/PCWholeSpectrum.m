function [pc,coeff]=PCWholeSpectrum(sps);
wspec=[];
% max_spec=max(max(sps.specs));
% thr=20*log10(max_spec)-80;
% sps.specs=20*log10(sps.specs);
% sps.specs=max(max(sps.specs,thr),thr);

[nstim,ndim]=size(sps.specs);
nfreq=sps.nfreq;
nws=sps.nms;

for i=1:nstim;
 u =reshape(sps.specs(i,:),nfreq,nws);
 u=u(:,1:10:end);
 [nf,nd]=size(u);
 wspec=[wspec;reshape(u,1,nf*nd)];
end;
u =reshape(sps.specs(end,:),nfreq,nws);
 u=u(:,1:10:end);
[pc,coeff]=princomp(wspec);

%u=zscore(wspec)*pc;
%u(1:10,1:10)
%coeff(1:10,1:10)