function HarmFun(x);
za=timefreq(x,44100,'stft');
ct=2;
n=length(za.t);
nslices=floor(length(za.t)/ct);
harms={}
slices=1;
hf=[]
ff=[];
dB=1e-1;
for start=1:ct/2:n;
    nct=ct;
    if start+ct>n;
        nct=n-start;
    end
    start+nct
    n
    size(za.spec)
    zs=za.spec(:,start:(start+nct));
    ec50=EC50(zs);
    if ec50<15;
        f=0;
        h=0;
        fm=0;
    else
       
       zzs=mean(zs,2);
       z=zzs-mean(zzs);
       nf=length(z);
       My=-1e6;
       mf=-1;
       [m,fm]=max(zzs);
       if fm==1;
           [m,fm]=max(zzs(3:end));
       end;
       d=[];
       for j=1:max(fm/5,1);
           fj=floor(fm/j)
           if zzs(fj)>zzs(fm)*dB;
               d=[d;fj];
           end;
       end;
       dm=min(d);
       
       mf=dm;
       f=mf;
       z=zzs-max(zzs)*dB;
       h=length(find(z(2*f:f:nf)>0));           
    end;
    harms{slices}.f=f;
    harms{slices}.h=h;
    
    slices=slices+1;
    
    hf=[hf ;start ec50 f];
    ff=[ff;start f fm];
    
    for k=1:h;
        hf=[hf;start ec50 f+k*f];
    end;
end;
    t=1:0.1:n;
 subplot(2,2,1)

 hhf=HarmFun2(x)
    
    plot(hf(:,1)/n*za.t(end),hf(:,3)/60*za.f(end),'+',t/n*za.t(end),interp1(ff(:,1),ff(:,2)/60*za.f(end),t,'cubic'),hhf(:,1)/n*za.t(end),hhf(:,2))
    
    axis([0 za.t(end) 0 za.f(end)])
    
    subplot(2,2,2)
   imagesc(za.t,za.f,za.spec);
   axis xy;
    
    
    subplot(2,2,3)

   
    plot(t,interp1(ff(:,1),ff(:,2)/60*za.f(end),t,'cubic'));
    subplot(2,2,4)
    t=1:0.1:slices;
    
    plot(ff(:,1),ff(:,3)/60*za.f(end))
