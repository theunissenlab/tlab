function ridge_mat=CreateChosenRidgeMatrix(psths,sps,nbin,tmin,tmax,deltat,pcs);

[nstim,npsth]=size(psths);
ridge_mat=[];
for psi=1:nstim;   
    psi
    psth=psths(psi,:);    
    spec=sps.specs(psi,:);
    
    stimSpectrum=reshape(spec,sps.nfreq,sps.nms);   
    psthRMat=ComputeSpectrumRidgeMatrix(psth,nbin,tmin,tmax,stimSpectrum,deltat,pcs);
    ridge_mat=[ridge_mat;psthRMat];
end;
