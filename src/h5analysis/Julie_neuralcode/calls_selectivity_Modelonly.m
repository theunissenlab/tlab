function [dsel, AvgDprime, StdDprime, Comptype, AvgRcallcat,StdRcallcat]=calls_selectivity_Modelonly(MatfilePath)
Cal = struct();
Res=load(MatfilePath);
StimType=unique(Res.VocType);
NstimType=length(StimType);
Nsections=length(Res.VocType);

%% Calculate mean and std spike rate for each call category
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

%% calculate index of selectivity based on average spike rate over call categories 
Ldsel=cumsum(NstimType);
Ldsel=Ldsel(end);
dsel=cell(Ldsel,1);
Comptype=cell(Ldsel,2);
dd=0;
for vt1=1:(NstimType-1)
    for vt2=(vt1+1):NstimType
        dd=dd+1;
        codevoc1=sum(double(char(StimType(vt1))));
        codevoc2=sum(double(char(StimType(vt2))));
        if codevoc1<=codevoc2
            dsel{dd}=(AvgRcallcat(vt1)-AvgRcallcat(vt2))/(AvgRcallcat(vt1)+AvgRcallcat(vt2));
            Comptype{dd,1}=[StimType(vt1) StimType(vt2)];
            Comptype{dd,2}= codevoc1*codevoc2;
        elseif codevoc1>codevoc2
            dsel{dd}=(AvgRcallcat(vt2)-AvgRcallcat(vt1))/(AvgRcallcat(vt1)+AvgRcallcat(vt2));
            Comptype{dd,1}=[StimType(vt2) StimType(vt1)];
            Comptype{dd,2}= codevoc1*codevoc2;
        end
    end
end
fprintf(1,'index of selectivity based on average spike rate over call categories: done\n');

%% Calculate all the within-dprimes
icomp = 1;
nfiles=length(Res.MeanRate);
Ldp=cumsum(nfiles);
Ldp=Ldp(end);
dprime=cell(Ldp,1);
AllComptype=cell(Ldp,1);
stim_mean_rate=cell2mat(Res.MeanRate);
spike_std_rate=cell2mat(Res.StdRate);
for nf1=1:nfiles-1
    for nf2=nf1+1:nfiles
        if (spike_std_rate(nf1) == 0 && spike_std_rate(nf2) == 0 )
            dprime{icomp} = 0;
       %work from here!!     AllComptype{icomp} = sum(double(char(Res.VocType(nf1)))) *sum(double(char(Res.VocType(nf2))));
        else
            codevoc1=sum(double(char(Res.VocType(nf1))));
            codevoc2=sum(double(char(Res.VocType(nf2))));
            if codevoc1<=codevoc2
                dprime{icomp} = 2*(stim_mean_rate(nf1)-stim_mean_rate(nf2))./sqrt(spike_std_rate(nf1)^2+spike_std_rate(nf2)^2);
                AllComptype{icomp} = codevoc1*codevoc2;
            elseif codevoc1>codevoc2
                dprime{icomp} = 2*(stim_mean_rate(nf2)-stim_mean_rate(nf1))./sqrt(spike_std_rate(nf1)^2+spike_std_rate(nf2)^2);
                AllComptype{icomp} = codevoc1*codevoc2;
            end
        end
        icomp = icomp +1;
    end
end
fprintf(1,'all the within-dprimes: done\n');

%% Calculate mean and std of dprime between call type
LComptype=length(cell2mat(Comptype(:,2)));
for vt=1:NstimType%this loop first complete the comparisons type list with the comparisons within the same group
    Comptype{LComptype+vt,1}=[StimType(vt) StimType(vt)];
    Comptype{LComptype+vt,2}= sum(double(char(StimType(vt))))*sum(double(char(StimType(vt))));
end

NComptype=length(cell2mat(Comptype(:,2)));
AvgDprime=cell(NComptype,1); 
StdDprime=cell(NComptype,1);
for ct = 1:NComptype
    cti=Comptype{ct,2};
    selector=(cell2mat(AllComptype)==cti);
    selectorInd=find(selector);
    AvgDprime{ct}=mean(cell2mat(dprime(selectorInd)));
    StdDprime{ct}=std(cell2mat(dprime(selectorInd)));
end
fprintf(1,'mean and std of dprime between call type: done\n');


%% Fit the neural responses to vocalization type and/or the spectro
fprintf(1, 'fitting the neural responses to vocalization type and/or the spectro\n');
[R2A, ModelPredict, LL, NEC, PValLRatio, h, NeuralResponse, voc, Best_nbPC]=Spectro_Neuro_model(MatfilePath);

%% Store Values
fprintf(1, 'Storing values\n');
Cal.subject=Res.subject;
Cal.dsel=dsel;
Cal.AvgDprime=AvgDprime;
Cal.StdDprime=StdDprime;
Cal.Comptype=Comptype;
Cal.AvgRcallcat=AvgRcallcat;
Cal.StdRcallcat=StdRcallcat;
Cal.Best_nbPC=Best_nbPC;
Cal.R2A=R2A;
Cal.ModelPredict = ModelPredict;
Cal.LogLikelihood = LL;
Cal.PValLRatio = PValLRatio;
Cal.SignifModelCompare = h;
Cal.NeuralResponse = NeuralResponse;
Cal.VocType=voc;


if ismac()
        [status username] = system('who am i');
        if strcmp(strtok(username), 'frederictheunissen')
            if strncmp('/auto/fdata/solveig',stim_name, 19)
            elseif strncmp('/auto/fdata/julie',stim_name, 17)
                calfilename = fullfile('/Users','frederictheunissen','Documents','Data','Julie','matfile',Res.subject,['Cal' Res.Site '.mat']);
            end
        elseif strcmp(strtok(username), 'elie')
            calfilename = fullfile('/Users','elie','Documents','MATLAB','data','matfile',Res.subject,['Cal' Res.Site '.mat']);
        end
else
    calfilename=fullfile('/auto','k6','julie','matfile',Res.subject,['Cal' Res.Site '.mat']);
end

save(calfilename, '-struct', 'Cal');
clear Res Cal  dsel AvgDprime StdDprime Comptype AvgRcallcat StdRcallcat R2A ModelPredict R2A LL NEC PValLRatio h NeuralResponse
%clear gamma_mutual_info gamma_noise_entropy totalentropy gamma_const
%rate_bandwidth rate_gamma fano_factor rate_info_biased rate_info_bcorr rate_info_stderr
end 
