function [pcs]=GetPCSpectrum(sps,deltat,overlap,nPC);

% Compute the pcs of the reduced sliding spectrum over deltat
% with overlap. 
% 
% INPUT : 
%        sps: the cell obtain using the function 
%        deltat: the size in ms of the sliding window 
%        overlap: in [0,1[ yielding the amount of overlap
%        nPC:the number of PC to return 
%  OUTPUT: 
%         pcs struct with 
%         .pc the actual pc matrix nparam*nPC matrix 
%         .nparam (=deltat*nfreq) the number of params total parsed
%         .mu : the overall mean 
%         .sigma : the overall std 
%     
%      

specs=sps.specs;
nms=sps.nms;
nfreq=sps.nfreq;
[nstims,nw]=size(specs);

if nw ~=nms*nfreq;
    disp('error: mismatch spectrums');
else

pcs={};
pcs.nparam=deltat*nfreq;
%pcs.nparam=nparam*nPC;

% we first create the overall matrix 
b_mat=[];
stim_v=[];
for stim=1:nstims;    
    spec=specs(stim,:);
    spec=reshape(spec,nfreq,nms);
    ov=1;
    while ov+deltat<nms;    
        sp=reshape(spec(:,ov:(ov+deltat-1)),1,deltat*nfreq);
        b_mat=[b_mat;sp];
        stim_v=[stim_v;stim ov-1];
        ov = ov + floor((1-overlap)*deltat);        
    end;
    rem=nms-ov;
    if rem>10;
        [nf,nr]=size(spec(:,ov:end));
        sp=reshape(spec(:,ov:end),1,nr*nfreq);
        sp=[sp,zeros(1,(deltat-nr)*nfreq)];        
        b_mat=[b_mat;sp];
        stim_v=[stim_v;stim ov-1];
    end;
end;
%pcs.b_mat=b_mat;
[z_mat,mu,sigma]=zscore(b_mat);
[pc,co]=princomp(b_mat);
pcs.pc=pc(:,1:nPC);
pcs.mu=mu;
pcs.sigma=sigma;
pcs.co=co(:,1:nPC);
pcs.stim_v=stim_v;
end;