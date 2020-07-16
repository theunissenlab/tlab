function [res,psi,psi2]=ParseUnit(unit,stim_dir,motograms_dir)
% general info
[path,bird,ext]=fileparts(unit.source_directory);
root_path='data';

site=unit.site;
depth=unit.rdepth;
electrode=unit.electrode;

% psth variable
tmin=0;% ms
tmax=600;% ms
nbin=20; % ms

deltat=100; % ms
overlap=0.5; % 50% overlap
nPC=20; % first npcs

nclasses=length(unit.classes);
for ncl=1:nclasses;
    % get the correct unit
    class=unit.classes{ncl};
    responses=unit.class_responses.(class);
    % get the psths and stims
    [psths,stims,stim_dir]=GetUnitPsth(responses,tmin,tmax,nbin);
    stim_dir = ['/Users/hsoula/Research/AudioNeuro/SolveigData/',stim_dir];
    stim_distance=[];
    stim_names={};    
    for i=1:length(stims)
        name=stims{i};
        ixd=find(name=='_');
        out=str2num(name((ixd+1):(end-1)));
        stim_distance=[stim_distance;out];
        stim_names{i}=name(1:(ixd-1));
    end;
    [nstims,npsth]=size(psths);
    [psths_n,mpsth,spsth]=NormalizePSTH(reshape(psths,nstims*npsth,1),npsth);
    % create directory
    warning('off','all');
    data_path=sprintf('%s/%s_%s',root_path,bird,site);
    mkdir(data_path);
    
    % compute spectrum
    spec_hash_name=sprintf('%s/specs_%s_%s_%d_%d_%s.mat',data_path,bird,site,depth,electrode,class);
    if ~exist(spec_hash_name,'file');
        sps=GetAllSpectrum(stims,stim_dir);
        max_spec=max(max(sps.specs));
        thr=20*log10(max_spec)-80;
        sps.specs=20*log10(sps.specs);
        sps.specs=max(max(sps.specs,thr),thr);
        
        save(spec_hash_name,'sps');
    else
        load(spec_hash_name);
    end;
    
    % get motograms 
    moto_hash_name=sprintf('%s/motos_%s_%s_%d_%d_%s.mat',data_path,bird,site,depth,electrode,class);
    if ~exist(moto_hash_name,'file');
        motos=GetMotograms(stims,stim_dir,motograms_dir);
        save(moto_hash_name,'motos');
    else
        load(moto_hash_name);
    end;
    
    %compute pcs
    pcs_hash_name=sprintf('%s/pcs_%s_%s_%d_%d_%s.mat',data_path,bird,site,depth,electrode,class);
    if ~exist(pcs_hash_name,'file');
        pcs=GetPCSpectrum(sps,deltat,overlap,nPC);
        save(pcs_hash_name,'pcs');
    else
        load(pcs_hash_name);
    end;
    
    
    ridge_hash_name=sprintf('%s/rdg_%s_%s_%d_%d_%s.mat',data_path,bird,site,depth,electrode,class);
    if ~exist(ridge_hash_name,'file');
        ridge_mat=CreateRidgeMatrix(psths,sps,nbin,tmin,tmax,deltat,pcs);
        save(ridge_hash_name,'ridge_mat');
    else
        load(ridge_hash_name);
    end;
    
    res={};
    
    
    
    [pc,co]=PCWholeSpectrum(sps);
    cop=co(:,1:20);
    y=mean(psths,2);
    res.mean_psths=y;
    res.dists=stim_distance;
    res.stim_names=stim_names;
    dists=[2,16,64,128,256];
    %    figure(1);
    [pred,b,r2,psi]=RidgePlot(y,cop,0.5,stim_distance,dists);
    res.all_spec_b=b;
    cop_2m=cop;
    for i=2:length(dists);
        stim_distance(i);
        fi=find(stim_distance==stim_distance(i));
        ys=y(fi);
        f2m=find(stim_distance==stim_distance(1));
        for i=1:length(fi);
            ij=fi(i);
            name=stim_names(ij);
            si=find(strcmp(stim_names,name));
            ni=intersect(f2m,si);
            cop_2m(ij,:)=cop(ni(1),:);            
        end
    end
  %  res.cop_2m=cop_2m;
 %   res.cop=cop;
    
     [pred_2m,b2m,r2_m,psi_m]=RidgePlot(y,cop_2m,0.5,stim_distance,dists);
    res.psi_m=psi_m;    
    
    %r2
    res.mean_prediction=r2;
    res.d2m_only=GetR2(psi{1}.ys,psi{1}.bsi(1)+psi{1}.ds*psi{1}.bsi(2:end));
    res.spec_psi=psi;
    cop=motos.co(:,1:20);
    y=mean(psths,2);
    dists=[2,16,64,128,256];
    
 %   figure(4);
    [pred,b,r2,psi2]=RidgePlot(y,cop,0.5,stim_distance,dists);
    %r2
    res.moto_mean_prediction=r2;
    res.d2m_moto_only=GetR2(psi2{1}.ys,psi2{1}.bsi(1)+psi2{1}.ds*psi2{1}.bsi(2:end));
    res.moto_psi=psi2;
    co=pcs.co;
    stim_v=pcs.stim_v;
    [nstim_v,nr]=size(stim_v);
    y=[];z=[];
    [nstim,tm]=size(psths);
    for i=1:nstim_v;
        sti=stim_v(i,1);
        out=stim_distance(sti);
        z=[z;out];
        ti=stim_v(i,2)+deltat;
        npst=floor(ti/nbin);
        npse=min(floor((ti+4*deltat)/nbin),tm);
        y=[y;mean(psths(sti,npst:npse))];
    end;
    
    
   % figure(2);
    [pred,b,r2]=RidgePlot(y,co,0.1,z,dists);
    %r2
    res.forward_prediction=r2;
    psths_d=zeros(nstim,30);
    
    for next=0:30;   
        z=[];id=[];
        nt=5;    
        nn=0;
        %next=10;
        y=[];
        for i=1:30;
            un=max(min(30,i+next),1);
            vn=max(min(30,i+nn),1);
            psths_d(:,i)=mean(psths(:,vn:un),2);
        end
        for i=1:nstim;
            out=stim_distance(i);
            z=[z; out*ones(30-nt,1)];
            y=[y;psths_d(i,(nt+1):end)'];
            id=[id;(i-1)*30+((nt+1):30)'];
        end;
      %  figure(3);
        [pred,b,r2]=RidgePlot(y,ridge_mat(id,1:20),0.5,z,dists);
        snext=sprintf('next_%d',next);
        res.(snext)=r2;
    end
    
    result_name=sprintf('%s/res_%s_%s_%d_%d_%s.mat',data_path,bird,site,depth,electrode,class);
    save(result_name,'res');    
end;