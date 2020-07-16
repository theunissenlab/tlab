function spl=fitWave(wavname,thr,freq_a)

w=wavread(wavname);

spec=htimefreq(w,44100,'stft');
splits=SplitSong(spec,thr);
spl={};
spl.w=w;
[path,wname,ext]=fileparts(wavname);
spl.wavname=wname;

spl.nw=0*w;
spl.r_moto={};
spl.splits={};
nl=length(splits);
for n=1:nl;
    wavf={};     
    [b_l,a_l]=FitSplit(splits{n},freq_a,n,nl);
    wavf.b=b_l;
    wavf.a=a_l;
    spl.r_moto{n}=wavf;
    spl.splits{n}=splits{n};
end;
