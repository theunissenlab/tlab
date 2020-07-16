function expl=GetAllData(bird,path);
expl={};
nsite=6;
nelectrode=16;

for site=1:nsite;
        dpath=sprintf('%s/%s_Site%d',path,bird,site);
        ssite=sprintf('s%d',site);
        expl{site}={};
    for electrode=1:nelectrode;
        dwild=sprintf('%s/res_%s_Site%d_*_%d_Mask*.mat',dpath,bird,site,electrode);
        fs=dir(dwild);
        nfs=length(fs);
        if nfs>1;
            fprintf('error: should be one file %s',dwild);
        elseif nfs==0;
                fprintf('no site for %d\n',site);            
        else
            %fs(1)
            u=expl{site};
            load(sprintf('%s/%s',dpath,fs(1).name));          
            %selectrode=sprintf('e%d',electrode);
            spec_p16dm_16mdata=GetR2(res.spec_psi{2}.ys,res.spec_psi{2}.bsi(1) + res.spec_psi{2}.ds*res.spec_psi{2}.bsi(2:end));            
            moto_p16dm_2mdata=GetR2(res.moto_psi{2}.ys,res.moto_psi{2}.bsi(1) + res.moto_psi{2}.ds*res.moto_psi{2}.bsi(2:end));
            
            spec_p16dm_2mdata=GetR2(res.psi_m{2}.ys,res.psi_m{2}.bsi(1) + res.psi_m{2}.ds*res.psi_m{2}.bsi(2:end));            
            %spec16my=res.spec_psi{2}.ys;
            %spec16md=res.spec_psi{2}.ds;
            
            expl{site}{electrode}=[res.d2m_only,res.d2m_moto_only,spec_p16dm_16mdata,moto_p16dm_2mdata,spec_p16dm_2mdata];
        end
    end;
end;
