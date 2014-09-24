function [fm,spec,res]=plotBGAspectrums(t,a,ta,b,tb)


y=smBGAs(t,24000,a,ta,b,tb);
tfrep=timefreq(y,44100,'stft');

tt=linspace(0,t,10);

%plot(interp1(ta,a,tt,'cubic'),interp1(tb,b,tt,'cubic'))
subplot(2,3,6)
plot(ta,a/max(a),tb,b/max(b))
subplot(2,3,1)
plot(y)
subplot(2,3,2)
plot(sum(tfrep.spec,2))
subplot(2,3,5)
plot(sum(tfrep.spec,1))
subplot(2,3,4)

imagesc(tfrep.t, tfrep.f, tfrep.spec);
axis xy;
v_axis = axis;
v_axis(1) = min(tfrep.t);
v_axis(2) = max(tfrep.t);
v_axis(3) = min(tfrep.f);
v_axis(4) = max(tfrep.f);
axis(v_axis);
xlabel('Time'), ylabel('Frequency');

spec=tfrep.spec;
subplot(2,3,3)
res=[];
for f=1:60;
    y=computeCOf(sum(spec,2),f);
    res=[res;f y];
    plot(res(:,1),res(:,2));
    drawnow;
end
[m,fm]=max(res(:,2));
fm=fm/60*8000;