function [psths_c,ridge_mat_c]=CreateChosenRidgeMatrix(chosen,psths,sps,nbin,tmin,tmax,deltat,pcs);

nchose=length(chosen);
psths_c=[];
ridge_mat_c=[];
for nc=1:nchose;
    psi=chosen(nc);
    psth=psths(psi,:);    
    psths_c=[psths_c;psth'];
    spec=sps.specs(psi,:);
    stimSpectrum=reshape(spec,sps.nfreq,sps.nms);
    psthRMat=ComputeSpectrumRidgeMatrix(psth,nbin,tmin,tmax,stimSpectrum,deltat,pcs);
    ridge_mat_c=[ridge_mat_c;psthRMat];
end;
