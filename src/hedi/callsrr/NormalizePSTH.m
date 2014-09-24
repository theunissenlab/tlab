function [npsth,mpsth,spsth]=NormalizePSTH(psth,nparam)

n=length(psth);
nstim=n/nparam;
m=[];
s=[];
for i=1:nparam;
    m=[m;mean(psth(i:nparam:end))];
    s=[s;std(psth(i:nparam:end))];
end;
me=mean(psth(1:nparam:end));
mpsth=m;
spsth=s;
npsth=0*psth;
for i=1:nstim;
   offset=1+(i-1)*nparam;
    npsth(offset:(offset+nparam-1)) = (psth(offset:(offset+nparam-1)) - me)./spsth;
end;