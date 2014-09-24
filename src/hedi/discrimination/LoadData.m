function [specs,motos,group]=LoadData(sex,wave_dir,NF_dir);

% [specs,motos,group]=LoadData(sex,wave_dir,NF_dir);

nms=GetLongestNF(sex,NF_dir);


%waves=[];
motos=[];
specs=[];
group=[];
nd=10;

for nbird=1:16;
    for ncall=1:16;
        birdname=sprintf('%cbird%dcall%d',sex,nbird,ncall);
        NF_name=sprintf('%s/NF2.%s.txt',NF_dir,birdname);
        wave_name=sprintf('%s/%s.wav',wave_dir,birdname);
        group=[group;nbird];
        moto=load(NF_name);
        [nti,nv]=size(moto);
        ntemp=zeros(nms,nv-1);
        ntemp(:,1)=-0.4;
        ntemp(:,2)=mean(moto(:,3));
        ntemp(:,3)=mean(moto(:,4));
        st=max(floor(nms/2-nti/2),1);
        en=st+nti-1;
        
        ntemp(st:en,:)=moto(1:end,2:end);        
        motos=[motos;reshape(ntemp,1,nms*(nv-1))];
        
        
        [out,samplerate]=wavread(wave_name);

        tfrep=timefreq(out,samplerate,'stft');
        ev=sum(tfrep.spec);
        tfrep.t=tfrep.t(ev>max(ev)/10);
        tfrep.spec=tfrep.spec(:,ev>max(ev)/10);
        tfrep.t=tfrep.t-tfrep.t(1);
        
        [nf,nt]=size(tfrep.spec);
        nstart=(nms-nt)/2;
        
        u=[zeros(nf,nstart),tfrep.spec,zeros(nf,nms-nstart-nt)];
        u=u(:,1:nd:end);
        [nf,ndec]=size(u);
        specs=[specs;reshape(u,1,nf*ndec)];
         
    end
end