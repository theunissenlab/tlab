function res=analyze(sex,specs,motos,groups);

wav_dir='../data';
NF_dir='../output';

%[specs,motos,groups]=LoadData(sex,wav_dir,NF_dir);

[zsp,sp_mu,sp_si]=zscore(specs);
[zmo,mo_mu,mo_si]=zscore(motos);


%[spec_pc,spec_co]=princomp(specs);
%[motos_pc,motos_co]=princomp(motos);
errorsSP=[];
errorsMO=[];

errorsSPPC=[];
errorsMOPC=[];


nDF=5;

meSP=[];
meMO=[];

meSPPC=[];
meMOPC=[];
res={};
for nPC=2:20;
    for ntrial=1:100;
        %    ntrial
        erMo=DFAs(zmo,groups,16,0.80,nPC,min(nPC,nDF));
        erSp=DFAs(zsp,groups,16,0.80,nPC,min(nPC,nDF));
        errorsSP=[errorsSP;erSp];
        errorsMO=[errorsMO;erMo];
        
        erMo=DFAs(zmo,groups,16,0.80,nPC,min(nPC,2));
        erSp=DFAs(zsp,groups,16,0.80,nPC,min(nPC,2));
        errorsSPPC=[errorsSPPC;erSp];
        errorsMOPC=[errorsMOPC;erMo];
        
    end;
    
    meSP=[meSP; nPC mean(errorsSP) std(errorsSP)/10];
    meMO=[meMO; nPC mean(errorsMO) std(errorsMO)/10];
    meSPPC=[meSPPC; nPC mean(errorsSPPC) std(errorsSPPC)/10];
    meMOPC=[meMOPC; nPC mean(errorsMOPC) std(errorsMOPC)/10];
    subplot(2,1,1)
    errorbar(meSP(:,1),meSP(:,2),meSP(:,3),'k-');
    hold on
    errorbar(meMO(:,1),meMO(:,2),meMO(:,3),'--r');
    drawnow
    subplot(2,1,2)   
    errorbar(meSPPC(:,1),meSPPC(:,2),meSPPC(:,3),'k-');
    hold on
    errorbar(meMOPC(:,1),meMOPC(:,2),meMOPC(:,3),'--r');
    drawnow
end

res.mespc5=meSP;
res.mesMO5=meMO;
res.mespc2=meSPPC;
res.mesMO2=meMOPC;
