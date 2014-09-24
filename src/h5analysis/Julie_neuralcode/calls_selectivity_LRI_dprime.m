function [calfilename]=calls_selectivity_LRI_dprime(MatfilePath)
%Whine vocalizations are not studied here except for average spike rate
%where the value is calculated for the units that have it.
Cal = struct();
Res=load(MatfilePath);
StimType=unique(Res.VocType);
RestrictAna=find(~strcmp(StimType, 'Wh'));% localize the stims we don't want to study here and get rid of them
if sum(strcmp(StimType, 'mlnoise'))%make sure that there are some mlnoise for that file. 
    RestrictAna(find(RestrictAna==find(strcmp(StimType, 'mlnoise'))))=[];
end
StimTypeR = StimType(RestrictAna);
NstimType=length(StimType);
NstimTypeR=length(StimTypeR);

%% Calculate mean and std spike rate (AvgRcallcat and StdRcallcat) for each call category
AvgRcallcat=zeros(NstimType,1);
StdRcallcat=zeros(NstimType,1);
for vt=1:NstimType
    st=StimType(vt);
    selector=strcmp(Res.VocType, st);
    selectorInd=find(selector);
    AvgRcallcat(vt)=mean(cell2mat(Res.MeanRate(selectorInd)));
    StdRcallcat(vt)=std(cell2mat(Res.MeanRate(selectorInd)));
end
fprintf(1,'mean and std spike rate for each call category: done\n');


%% Calculate SSI Semantic Selectivity Index and LogRatio Index LRI
AvgRcallother=zeros(NstimTypeR,1);
AvgRcallself=zeros(NstimTypeR,1);
for vt=1:NstimTypeR
    selectorI=RestrictAna(vt);
    AvgRcallself(vt)=AvgRcallcat(selectorI);
    selectorOther=RestrictAna;
    selectorOther(vt)=[];
    AvgRcallother(vt)=nanmean(AvgRcallcat(selectorOther));
end
SSIcat=(AvgRcallself-AvgRcallother)./(AvgRcallself + AvgRcallother);
AvgRcallself(find(AvgRcallself==0))=min(AvgRcallself(find(isfinite(log2(AvgRcallself)))))/100;% to get rid of infinite values of LRI I choose to replace nule values of spike rate by 1/100 of the min value of spike rate for that cell
LRIcat=log2(AvgRcallself./AvgRcallother);

%% Calculate dprimes (d-prime) between every single section of vocalization and section of mlnoise
if sum(strcmp(StimType, 'mlnoise'))
    icomp = 0;
    selector_mln=find(strcmp(Res.VocType, 'mlnoise'));
    selector_other=find(~strcmp(Res.VocType, 'Wh'));%suppresses indices of selector_voc that corresponds to Whines 'Wh' vocalizations
    selector_other=setdiff(selector_other,selector_mln);
    nOther=length(selector_other);
    nmln=length(selector_mln);
    Ldpmln=cumsum(1:(nmln-1));
    Ldpmln=Ldpmln(end);
    Ldp=nmln*nOther+Ldpmln;
    dprime=cell(Ldp,1);
    AllComptype=cell(Ldp,1);
    stim_mean_rate=cell2mat(Res.MeanRate);
    spike_std_rate=cell2mat(Res.StdRate);
    for nf1=1:nOther
        nO=selector_other(nf1);%% problem when mlnoise because he calculate twice and 
        for nf2=1:nmln
            nM=selector_mln(nf2);
            icomp = icomp +1;
            if (spike_std_rate(nO) == 0 && spike_std_rate(nM) == 0 )
                dprime{icomp} = 0;
                AllComptype{icomp} = Res.VocType{nO};
            else
                dprime{icomp} = 2*(stim_mean_rate(nO)-stim_mean_rate(nM))./sqrt(spike_std_rate(nO)^2+spike_std_rate(nM)^2);
                AllComptype{icomp} = Res.VocType{nO};
            end
        end
    end
    for n1=1:nmln-1
        nn1=selector_mln(n1);
        for n2=n1+1:nmln
            nn2=selector_mln(n2);
            icomp=icomp+1;
            if (spike_std_rate(nn1) == 0 && spike_std_rate(nn2) == 0 )
                dprime{icomp} = 0;
                AllComptype{icomp} = Res.VocType{nn1};
            else
                dprime{icomp} = 2*(stim_mean_rate(nn1)-stim_mean_rate(nn2))./sqrt(spike_std_rate(nn1)^2+spike_std_rate(nn2)^2);
                AllComptype{icomp} = Res.VocType{nn1};
            end
        end
    end

    fprintf(1,'all the dprimes between every single stim and mln: done\n');
else
    dprime=cell(0);
    AllComptype=cell(0);
    fprintf(1,'No dprimes between every single stim and mln because there are no mln for that file\n');
end

%% Calculate mean, std and CV of dprime
ComptypeU=unique(AllComptype);
LComptype=length(ComptypeU);
MeanDprime=zeros(LComptype,1);
STDDprime=zeros(LComptype,1);
dprimeM=cell2mat(dprime);
for vt=1:LComptype
    MeanDprime(vt)=nanmean(dprimeM(find(strcmp(AllComptype,ComptypeU(vt)))));
    STDDprime(vt)=nanstd(dprimeM(find(strcmp(AllComptype,ComptypeU(vt)))));
end
    
fprintf(1,'mean and std of dprime between call type: done\n');



%% Store Values
fprintf(1, 'Storing values\n');
Cal.subject=Res.subject;
Cal.StimType=StimType;
Cal.StimTypeR=StimTypeR;
Cal.AvgDprime=MeanDprime;
Cal.StdDprime=STDDprime;
Cal.AllComptype=AllComptype;
Cal.AvgRcallcat=AvgRcallcat;
Cal.StdRcallcat=StdRcallcat;
Cal.SSIcat=SSIcat;
Cal.LRIcat=LRIcat;
Cal.Dprime=dprime;

if ismac()
        [status username] = system('who am i');
        if strcmp(strtok(username), 'frederictheunissen')
            if strncmp('/auto/fdata/solveig',stim_name, 19)
            elseif strncmp('/auto/fdata/julie',stim_name, 17)
                calfilename = fullfile('/Users','frederictheunissen','Documents','Data','Julie','matfile',Res.subject,['LRI_DP' Res.Site '.mat']);
            end
        elseif strcmp(strtok(username), 'elie')
            calfilename = fullfile('/Users','elie','Documents','MATLAB','data','matfile',Res.subject,['LRI_DP' Res.Site '.mat']);
        end
else
    calfilename=fullfile('/auto','k6','julie','matfile',Res.subject,['LRI_DP' Res.Site '.mat']);
end

save(calfilename, '-struct', 'Cal');
clear Res dprime Cal

end 
