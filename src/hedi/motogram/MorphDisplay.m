function MorphDisplay(mgStart,mgEnd,nstep,out);

md=MgramMorph(mgStart,mgEnd,nstep);
mt=MgramMorphTimber(mgStart,mgEnd,nstep);
ms=MgramMorphStyle(mgStart,mgEnd,nstep);
MorphWave(md,sprintf('di-%s',out));
MorphWave(mt,sprintf('ti-%s',out));
MorphWave(ms,sprintf('si-%s',out));
ind=[1,2:4:(nstep-1),nstep];
nind=length(ind)
n=1;
for i=ind
    d=MgramWave(md{i});
    sd=timefreq(d,44100,'stft');   
    t=MgramWave(mt{i});
    st=timefreq(t,44100,'stft');   
    
    s=MgramWave(ms{i});
    ss=timefreq(s,44100,'stft');   
    
    subplot(3,nind,n);
    imagesc(st.t,st.f,st.spec);axis xy
    
    subplot(3,nind,n+nind);
    imagesc(sd.t,sd.f,sd.spec);axis xy
    
    subplot(3,nind,n+2*nind);
    imagesc(ss.t,ss.f,ss.spec);axis xy
    
    n=n+1;
end