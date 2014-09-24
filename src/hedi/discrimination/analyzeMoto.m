function resAll=analyzeMoto(moto,group,nDF);

resAll={};
resAll.tim=[];
resAll.moto=[];
resAll.syr=[];

[ncall,nd]=size(moto);
[zmo,zmu,zsi]=zscore(moto);
rank=floor(nd/4);
mtim=moto;
msyr=moto;

syr_idx=[1:rank]
tim_idx=[(rank+1):(3*rank)]
msyr(:,syr_idx)=repmat(zmu(syr_idx),ncall,1);

mtim(:,tim_idx)=repmat(zmu(tim_idx),ncall,1);


for nPC=nDF:30;   
   errMoto=DFAalls(zscore(moto),group,16,nPC,nDF);
   errTim=DFAalls(zscore(mtim),group,16,nPC,nDF);
   errSyr=DFAalls(zscore(msyr),group,16,nPC,nDF);
   
   resAll.tim=[resAll.tim; nPC errTim];
   resAll.moto=[resAll.moto; nPC errMoto];
   resAll.syr=[resAll.syr; nPC errSyr];
   plot(resAll.moto(:,1),resAll.moto(:,2),resAll.tim(:,1),resAll.tim(:,2),resAll.syr(:,1),resAll.syr(:,2));
   drawnow;
end

   