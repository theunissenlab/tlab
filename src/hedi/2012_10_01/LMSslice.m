function y=LMSslice(start,target_slice,tmax,nslices,ns);

c=start(1);
p=start(2);
k=start(3)*1e7;
b=start(4)*1e4;
x0=start(5);
y0=start(6);

nfilter=10;
A=start(7:16);
fsample=44100;
f=linspace(0,6000,nfilter);
F=[f 8000 10000  fsample/2];
A=[A,0.00,0.0, 0.0];
ff=interp1(F,A,0:10:(fsample/2));
ZA=invfreqz(ff,(0:10:(fsample/2))/fsample*2*pi,100,1);
%ZA=invfreqz(A,linspace(0,2*fmax/fsample*pi,length(A)),10,1);

out=sms(tmax,c,p,k,b,x0,y0);


%out2=out;
out2=filter(ZA,[1],out);
out2=filter(ZA,[1],out2);

spec=timefreq(out2,44100,'stft');
final=spec.spec;
%final=final/max(final);
v=sum(final,2);
w=sum(target_slice,2);
v=v/max(v);
w=w/max(w);
y=sum((v-w).^2)

 imagesc(final);
 drawnow