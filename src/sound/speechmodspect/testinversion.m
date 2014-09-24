% Test Code..... as of December 1, 1995.
len=2000;
window = [sin((1:200)/200*pi/2) ...
               ones(1,len-200-200) ...
               cos(-(1:200)/200*pi/2)];
input = sin((1:2000)*2*pi/10) .* window;
spec=abs(ComplexSpectrum(input,64,256));

figure;
imagesc(spec);

y=SpectrumInversion(spec,64,256,100);

figure;
plot(input);
hold on;
plot(y,'k');

specy = abs(ComplexSpectrum(input,64,256));
figure;
imagesc(specy);


