function re=RMSA(X);

a=X(1:10);
b=0.0105*ones(1,10);%BB;%X(11:20);
wavespectarget=wavread('01M03A_281106_012.wav');
tfrep=timefreq(wavespectarget,44100,'stft');
t=tfrep.t(end);

nt=10;
tf=linspace(0,t,nt);

nfilter=10;

fmax=7000; fsample=44100;
f=linspace(0,7000,nfilter);
F=[f 8000 10000  fsample/2];

A=[0 0 0 5000 10000 50000 30000 5000 0 0 0.01 0.0 0.0];

ff=interp1(F,A,0:10:(fsample/2));
ZA=invfreqz(ff,(0:10:(fsample/2))/fsample*2*pi,100,1);

y=smBGAs(t,24000,a',tf',b',tf');

y=filter(ZA,1,y);
y=filter(ZA,1,y);
tfrepy=timefreq(y,44100,'stft');
spe1=tfrep.spec(:); spe2=tfrepy.spec(:);
rho=mean((spe1-mean(spe1)).*(spe2-mean(spe2))/(std(spe1)*std(spe2)))
re=1-rho^2;
subplot(1,2,1)
imagesc(tfrepy.t, tfrepy.f, tfrepy.spec);
axis xy;
v_axis = axis;
v_axis(1) = min(tfrep.t);
v_axis(2) = max(tfrep.t);
v_axis(3) = min(tfrep.f);
v_axis(4) = max(tfrep.f);
axis(v_axis);
xlabel('Time'), ylabel('Frequency');

subplot(1,2,2)
imagesc(tfrep.t,tfrep.f, tfrep.spec);
axis xy;
v_axis = axis;
v_axis(1) = min(tfrep.t);
v_axis(2) = max(tfrep.t);
v_axis(3) = min(tfrep.f);
v_axis(4) = max(tfrep.f);
axis(v_axis);
xlabel('Time'), ylabel('Frequency');

drawnow