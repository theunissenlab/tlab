function y=plotBGAspectrums(t,a,ta,b,tb)


y=smBGAs(t,24000,a,ta,b,tb);

tfrep=timefreq(y,44100,'stft');
imagesc(tfrep.t, tfrep.f, tfrep.spec);
axis xy;
v_axis = axis;
v_axis(1) = min(tfrep.t);
v_axis(2) = max(tfrep.t);
v_axis(3) = min(tfrep.f);
v_axis(4) = max(tfrep.f);
axis(v_axis);
xlabel('Time'), ylabel('Frequency');