function resAll=analyzeAll(moto,sp,group,nDF);

resAll={};
resAll.sp=[];
resAll.moto=[];

for nPC=nDF:30;
   errSP=DFAalls(zscore(sp),group,16,nPC,nDF);
   errMoto=DFAalls(zscore(moto),group,16,nPC,nDF);
   resAll.sp=[resAll.sp; nPC errSP];
   resAll.moto=[resAll.moto; nPC errMoto];
   plot(resAll.sp(:,1),resAll.sp(:,2),'k-',resAll.moto(:,1),resAll.moto(:,2),'--r');
   drawnow;
end

   