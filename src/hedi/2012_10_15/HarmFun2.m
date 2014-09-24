function hf=HarmFun2(x);
za=timefreq(x,44100,'stft');
ct=5;
n=length(za.t);
nslices=floor(length(za.t)/ct);
harms={}
slices=1;
hf=[]
ff=[];
dB=1e-1;
for start=1:ct:n;
    nct=ct;
    if start+ct>n;
        nct=n-start;
    end;
    
    zs=za.spec(:,start:(start+nct));
    ec50=EC50(zs);
    if ec50<15;
        f=0;
        h=0;
        fm=0;
    else
        t=za.t(start:(start+nct));
        fr=za.f;
        f=fundamental_estimator(t,fr,zs);
        f=mean(f);               
    end
    hf=[hf;start f];
end;
