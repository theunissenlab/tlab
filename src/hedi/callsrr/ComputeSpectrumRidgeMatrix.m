function psthRMat=ComputeSpectrumRidgeMatrix(psth,nbin,fmin,fmax,stimSpectrum,deltat,pcs);

% Compute the ridge matrix for one psth and one stimulus 
% for each value of the psth put the -delta_t in time of the spectrum
% INPUT:
%         psth: the psth computed with fmin,fmax,nbin sent also in parameters
%         
%         stimSpectrum: the spectrum of the stim that start at fmin 
%         
%         deltat: the time windonw in ms (same as the pcs)
%         
%         pcs: matrix nparam by npc ; nparam/deltat should yield the frequency bands 
%         
%  OUTPUT:
%         psthRMat: the ntimebin by nPC matrix corresponding to the ridge regression matrix 
    

ntimebin=length(psth);
psthRMat=[];
pc=pcs.pc;
mu=pcs.mu;

%sigma=pcs.sigma;
nmax=100;
%me=max(max(stimSpectrum));
%sstimSpectrum=stimSpectrum/me;
%mu=0;
sigma=1;
[nfreq,nms]=size(stimSpectrum);
i=1;
u=1;
for t=-fmin:nbin:fmax;
    if i<=ntimebin;
        if t<0;
            w0=(mu)./sigma*pc;            
            psthRMat=[psthRMat;w0];
            u=u+1;
        elseif t<deltat;
            sp=stimSpectrum(:,1:t);       
            spc=reshape(sp,1,t*nfreq);
            nc=(deltat-t)*nfreq;
            rw=[mu(1:nc),spc];            
            w1=((rw-mu)./sigma)*pc;          
            psthRMat=[psthRMat;w1];
            u=u+1;
        elseif t<nms;
   %         t
            sp=stimSpectrum(:,(floor(t-deltat)+1):t);
            rsp=reshape(sp,1,deltat*nfreq);             
            w2=((rsp-mu)./sigma)*pc;            
            psthRMat=[psthRMat;w2];
            u=u+1;
        elseif t<nms+deltat;
            sp=stimSpectrum(:,(t-deltat+1):end);             
            spc=reshape(sp,1,(nms-t+deltat)*nfreq);
            nc=(nms-t+deltat)*nfreq+1;
            nmu=mu(nc:end);
            rw=[spc,nmu];
            w3=((rw-mu)./sigma)*pc;
            psthRMat=[psthRMat;w3];
            u=u+1;
        else
%            fprintf('T %f\n',t);
            w4=((mu)./sigma)*pc;
            psthRMat=[psthRMat;w4];
            u=u+1;
        end;
    end;
    i=i+1;     
end;
% for t=-fmin:nbin:0; % only zeros to be send here 
%     w0=((zeros(1,deltat*nfreq)-mu)./sigma)*pc(:,:);
%     psthRMat=[psthRMat;w0];
% end;

% for t=1:nbin:(deltat-1);
%     sp=stimSpectrum(1:t,:);
%     spc=reshape(sp,1,t*nfreq);
%     rw=[zeros(1,(deltat-t)*nfreq),spc];
%     w1=((rw-mu)./sigma)*pc(:,:);
%     psthRMat=[psthRMat;w1];
% end;
% 
% for t=deltat:nbin:nms;
%     sp=stimSpectrum((floor(t-deltat)+1):t,:);
%     rsp=reshape(sp,1,deltat*nfreq);
%     w=((rsp-mu)./sigma)*pc(:,:);
%     psthRMat=[psthRMat;w];
% end;
% 
% for t=(nms+1):nbin:(nms+deltat);
%      sp=stimSpectrum((t-deltat+1):end,:);        
%      spc=reshape(sp,1,(nms-t+deltat)*nfreq);
%      rw=[spc,zeros(1,(deltat+t-nms-deltat)*nfreq)];
%      w3=((rw-mu)./sigma)*pc(:,:);
%      psthRMat=[psthRMat;w3];
% end
% 
% for t=(nms+deltat+nbin):nbin:fmax;
%     w4=((zeros(1,deltat*nfreq)-mu)./sigma)*pc(:,:);
%     psthRMat=[psthRMat;w4];
% end;