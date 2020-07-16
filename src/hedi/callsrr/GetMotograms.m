function motos=GetMotograms(stims,stim_dir,motograms_dir);
% Return the spectrums of all stims as output by timefreq
% INPUT:
%       stims: the vector of all stims used nstims=length(stims)
%       stim_dir: where the wave files are 
%       motograms_dir : where the motograms are 
%       
% OUTPUT:
%       motos:
%          .motos: a nstims*4*ntime matrix 
%          .nms: length of the longest call in ms 
%           


% first compute the longest call
nms=GetLongest(stims,stim_dir);
nstims=length(stims);
motos={};
motos.nms=nms;
motos.m=[];
motos.name={};
i=1;
for ns=1:nstims;
    % load wave    
    name=stims{ns};
    ixd=find(name=='_');
    dist=str2num(name((ixd+1):(end-1)));
    if dist==2;
        call=name(1:(ixd-1));
        motos.name{i}=call;
        i=i+1;
        new_name=sprintf('%cbird%s',call(1),call(2:end));
        motoname=sprintf('%s/NF2.%s.txt',motograms_dir,new_name);
        mo=load(motoname);
        [nti,nv]=size(mo);
        ntemp=zeros(nms,nv-1);
        ntemp(:,1)=-0.4;
        ntemp(:,2)=mean(mo(:,3));
        ntemp(:,3)=mean(mo(:,4));
        st=floor(nms/2-nti/2);
        en=st+nti-1;
        ntemp(st:en,:)=mo(:,2:end);
        motos.m=[motos.m;reshape(ntemp,1,nms*(nv-1))];
    end;
end;

[z,mu,sigma]=zscore(motos.m);
motos.mu=mu;
motos.sigma=sigma;
[motos.pc,co]=princomp(motos.m);
motos.co=[];
for ns=1:nstims;
    % load wave    
    name=stims{ns};
    ixd=find(name=='_');
    dist=str2num(name((ixd+1):(end-1)));
    call=name(1:(ixd-1));
    for j=1:length(motos.name);
        %motos.name(j)
        %co(j,1)
        if strcmp(motos.name(j),call)
            motos.co=[motos.co;co(j,:)];
        end
    end
end

        