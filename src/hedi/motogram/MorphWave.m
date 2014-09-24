function s=MorphWave(morphs,pre_fix)
dur=floor(200*44100/1000);
n=length(morphs);
s=zeros(1,dur);
for i=1:n;
    s=[s,MgramWave(morphs{i}),zeros(1,dur)];
end;

wavwrite(s,44100,sprintf('%s.wav',pre_fix));