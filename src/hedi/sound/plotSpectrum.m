function [tfrep]=plotSpectrum(fname);

[out,samplerate,nbits]=wavread(fname);

tfrep=timefreq(out,samplerate,'stft');
ev=sum(tfrep.spec);
tfrep.t=tfrep.t(ev>max(ev)/10);
tfrep.spec=tfrep.spec(:,ev>max(ev)/10);

imagesc(tfrep.t, tfrep.f, tfrep.spec);
axis xy;
v_axis = axis;
v_axis(1) = min(tfrep.t);
v_axis(2) = max(tfrep.t);
v_axis(3) = min(tfrep.f);
v_axis(4) = max(tfrep.f);
axis(v_axis);
xlabel('Time'), ylabel('Frequency');
