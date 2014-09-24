function [dsel, AvgDprime, StdDprime, Comptype, AvgRcallcat,StdRcallcat]=calls_selectivity(MatfilePath)
confusioncallcategory=2; %set this value to 1 if you want confusion matrix calculated based on call categories and not per stim 
Plotfig=0; %set this to 1 if you want to see confusion matrices 
Cal = struct();
Res=load(MatfilePath);
StimType=unique(Res.VocType);
NstimType=length(StimType);
nfiles=length(Res.VocType);

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

%% calculate index of selectivity (dsel) based on average spike rate calculated just above over call categories 
Ldsel=cumsum(1:(NstimType-1));
Ldsel=Ldsel(end);
dsel=zeros(Ldsel,1);
Comptype=cell(Ldsel,2);
dd=0;
for vt1=1:(NstimType-1)
    for vt2=(vt1+1):NstimType
        dd=dd+1;
        codevoc1=sum(double(char(StimType(vt1))));
        codevoc2=sum(double(char(StimType(vt2))));
        if codevoc1<=codevoc2
            dsel(dd)=(AvgRcallcat(vt1)-AvgRcallcat(vt2))/(AvgRcallcat(vt1)+AvgRcallcat(vt2));
            Comptype{dd,1}=[StimType(vt1) StimType(vt2)];
            Comptype{dd,2}= codevoc1*codevoc2;
        elseif codevoc1>codevoc2
            dsel(dd)=(AvgRcallcat(vt2)-AvgRcallcat(vt1))/(AvgRcallcat(vt1)+AvgRcallcat(vt2));
            Comptype{dd,1}=[StimType(vt2) StimType(vt1)];
            Comptype{dd,2}= codevoc1*codevoc2;
        end
    end
end
fprintf(1,'index of selectivity based on average spike rate over call categories: done\n');
Dsel=struct();
Dsel.values=dsel;
Dsel.comptype=Comptype;

%% Calculate SSI Semantic Selectivity Index and LogRatio Index LRI
AvgRcallother=zeros(NstimType,1);
AvgRcallI=zeros(NstimType,1);
for vt=1:NstimType
    st=StimType(vt);
    selector=strcmp(Res.VocType, st);
    selectorInd=find(selector);
    selectorIndOther=find(~selector);
    AvgRcallI(vt)=mean(cell2mat(Res.MeanRate(selectorInd)));
    AvgRcallother(vt)=mean(cell2mat(Res.MeanRate(selectorIndOther)));
end
SSIcat=(AvgRcallI-AvgRcallother)./(AvgRcallI + AvgRcallother);
AvgRcallI(find(AvgRcallI==0))=min(AvgRcallI(find(isfinite(log2(AvgRcallI)))))/10;
LRIcat=log2(AvgRcallI./AvgRcallother);

%% Calculate dprimes (d-prime) between every single sections
icomp = 1;
Ldp=cumsum(1:(nfiles-1));
Ldp=Ldp(end);
dprime=cell(Ldp,1);
AllComptype=zeros(Ldp,2);
stim_mean_rate=cell2mat(Res.MeanRate);
spike_std_rate=cell2mat(Res.StdRate);
for nf1=1:nfiles-1
    for nf2=nf1+1:nfiles
        codevoc1=sum(double(char(Res.VocType(nf1))));
        codevoc2=sum(double(char(Res.VocType(nf2))));
        if (spike_std_rate(nf1) == 0 && spike_std_rate(nf2) == 0 )
            dprime{icomp} = 0;
            AllComptype(icomp,1) = codevoc1*codevoc2;
            if codevoc1==codevoc2
                AllComptype(icomp,2)=1;
            end
        elseif codevoc1<codevoc2
            dprime{icomp} = 2*(stim_mean_rate(nf1)-stim_mean_rate(nf2))./sqrt(spike_std_rate(nf1)^2+spike_std_rate(nf2)^2);
            AllComptype(icomp,1) = codevoc1*codevoc2;
        elseif codevoc1==codevoc2
            dprime{icomp} = 2*(stim_mean_rate(nf1)-stim_mean_rate(nf2))./sqrt(spike_std_rate(nf1)^2+spike_std_rate(nf2)^2);
            AllComptype(icomp,1) = codevoc1*codevoc2;
            AllComptype(icomp,2) = 1;
        elseif codevoc1>codevoc2
            dprime{icomp} = 2*(stim_mean_rate(nf2)-stim_mean_rate(nf1))./sqrt(spike_std_rate(nf1)^2+spike_std_rate(nf2)^2);
            AllComptype(icomp,1) = codevoc1*codevoc2;
        end
        icomp = icomp +1;
    end
end
fprintf(1,'all the dprimes between every single section: done\n');

%% Calculate mean, std and CV of dprime between/within call type
LComptype=length(cell2mat(Comptype(:,2)));
for vt=1:NstimType%this loop first complete the comparisons type list with the comparisons within the same group
    Comptype{LComptype+vt,1}=[StimType(vt) StimType(vt)];
    Comptype{LComptype+vt,2}= sum(double(char(StimType(vt))))*sum(double(char(StimType(vt))));
end

%select the within-category dprime values and between categories dprime
%values
wi_cat_selector=find(AllComptype(:,2));
be_cat_selector=find(~AllComptype(:,2));
Mean_DP_wi=mean(abs(cell2mat(dprime(wi_cat_selector))));
Std_DP_wi=std(abs(cell2mat(dprime(wi_cat_selector))));
Mean_DP_be=mean(abs(cell2mat(dprime(be_cat_selector))));
Std_DP_be=std(abs(cell2mat(dprime(be_cat_selector))));
CV_DP_wi=Std_DP_wi/Mean_DP_wi;
CV_DP_be=Std_DP_be/Mean_DP_be;

NComptype=length(cell2mat(Comptype(:,2)));
AvgDprime=cell(NComptype,1); 
StdDprime=cell(NComptype,1);
for ct = 1:NComptype
    cti=Comptype{ct,2};
    selector=(AllComptype(:,1)==cti);
    selectorInd=find(selector);
    AvgDprime{ct}=mean(abs(cell2mat(dprime(selectorInd))));
    StdDprime{ct}=std(abs(cell2mat(dprime(selectorInd))));
end
fprintf(1,'mean and std of dprime between call type: done\n');

%% Create the data set for the type of comparison you want:

if confusioncallcategory==1
    % here with Selector = type of stim (AggC, DisC Tet...)
    SpikeTrains=cell(NstimType,1);
    %Sections_len=zeros(Nsections*length(Res.Trials{1},1));
    Sections_len=zeros(NstimType,1);

    for vt=1:NstimType
        st=StimType(vt);
        selector=strcmp(Res.VocType, st);
        selectorIndMC=find(selector);
        NselectorInd=length(selectorIndMC);
        %Sections_len{vt}=zeros(NselectorInd*length(Res.Trials{selectorInd(1)}));
        ISec = 1;
        for NSI = 1:NselectorInd
            SI=selectorIndMC(NSI);
            SpikeTrains{vt} = [SpikeTrains{vt} ; Res.Trials{SI}];
            NSec = length(Res.Trials{SI});
            Sections_len(vt)=min(Res.SectionLength(SI));
            ISec=ISec+NSec;
        end
    end
elseif confusioncallcategory==0    
%% Just run the regular confusion Matrix code for each stim
SpikeTrains=Res.Trials;
Sections_len = Res.SectionLength;%here stim_len is in ms
nfiles=length(Res.VocType);
selectorIndMC=1:nfiles;
elseif confusioncallcategory==2
    %% select the first cut of each vocalization
    %% Organize the files per call category to obtain an ordered confusion matrix
    selector1=find(Res.Cut_orders==1);
    %selector2=find(strcmp(Res.VocType, 'Ag'));
    %selector3=find(strcmp(Res.VocType, 'DC'));
    %selectorIndMC=[intersect(selector1, selector2); intersect(selector1, selector3)];
    selectorIndMC = 1:length(Res.Cut_orders);
    SpikeTrains_temp=Res.Trials(selectorIndMC);
    SpikeTrains=cell(size(SpikeTrains_temp));
    Sections_len_temp = Res.SectionLength(selectorIndMC);
    Sections_len = zeros(size(Sections_len_temp));
    Res_Indices = zeros(size(selectorIndMC));
    nfiles = length(Sections_len_temp);
    cc=0;
    VocTypeList=Res.VocType(selectorIndMC);
    VocTypeSel=cell(size(VocTypeList));
        for vt=1:NstimType
            st=StimType(vt);
            selector=strcmp(VocTypeList, st);
            selectorInd=find(selector);
            NselectorInd=length(selectorInd);
            for cc_temp=1:NselectorInd
                cc=cc+1;
                SpikeTrains{cc}=SpikeTrains_temp{selectorInd(cc_temp)};
                Sections_len(cc)=Sections_len_temp(selectorInd(cc_temp));
                VocTypeSel{cc}=VocTypeList{selectorInd(cc_temp)};
                Res_Indices(cc)=selectorIndMC(selectorInd(cc_temp));
            end
        end
end

%% Set up parameters for all the Info calculations    
% Window sizes for calculations that depend on window length
winSize = [2 3 4 200];    % Window size for static Info calculation and for ideal observer
% winSize = [50];    % Window size for static Info calculation and for ideal observer
ns = length(winSize);

% Initialize all return values to zero
%***Ideal Observer Calculations (by Trial Template)***
percorrectB = zeros(1,ns);
mi_confusionB = zeros(1,ns);

percorrectCT = zeros(1,ns);
mi_confusionCT = zeros(1,ns);

zdistTB = zeros(1,ns);              % Uses std for other for both self and other in KL Divergence calc
pzdistTB = zeros(1,ns);
mizdistTB = zeros(1,ns);
mizdistncTB = zeros(1,ns);

zdistSB = zeros(1,ns);              % Recenters "other" distrubtions to zero to minimize additive variance across songs
pzdistSB = zeros(1,ns);             % and collapses distances across all trials of a single stimulus
mizdistSB = zeros(1,ns);
mizdistncSB = zeros(1,ns);
confusionMatrix = cell(1,ns);


%***Gamma Model Calculations***
gamma_mutual_info = 0;
gamma_noise_entropy = 0;
gamma_spike_rate = 0;
gamma_const = 0;
rate_bandwidth = 0;
rate_gamma = 0;
fano_factor = 0;

%***Rate Information Calculations***
rate_info_biased = zeros(1,ns);
rate_info_bcorr = zeros(1,ns);
rate_info_stderr = zeros(ns,2);

%% Calculations!!
if nfiles
    % With VR_distanceB
    for is=1:ns
        fprintf(1, 'Calculating confusion matrix with winSize = %d\n', winSize(is));
        % Reinsert this if using the confusion matrix in version B
        if confusioncallcategory==1
            [pc, mi_conf, zdT, pzdT, mi_zdT, mi_zdT_nc, zdS, pzdS, mi_zdS, mi_zdS_nc, confusion_matrix] = info_distanceB(NstimType, SpikeTrains , Sections_len, @VR_distanceB, winSize(is),Plotfig);
        else
            [pc, mi_conf, zdT, pzdT, mi_zdT, mi_zdT_nc, zdS, pzdS, mi_zdS, mi_zdS_nc, confusion_matrix] = info_distanceB(nfiles, SpikeTrains , Sections_len, @VR_distanceB, winSize(is), Plotfig);
        end

        % creating the compiled confusion matrix per call category
        confusion_matrix_CallType = zeros(NstimType, NstimType);
        for vtR=1:NstimType
            stR=StimType(vtR);
            selectorR=strcmp(VocTypeSel, stR);
            selectorIndR=find(selectorR);
            for vtC = 1:NstimType
                stC=StimType(vtC);
                selectorC=strcmp(VocTypeSel, stC);
                selectorIndC=find(selectorC);
                confusion_matrix_CallType(vtR,vtC)=sum(sum(confusion_matrix(selectorIndR, selectorIndC)));
            end
        end
        
        mi_confCT=info_matrix(confusion_matrix_CallType);
        pcCT=sum(diag(confusion_matrix_CallType));
        if Plotfig==1
            figure(5);
            imagesc(confusion_matrix_CallType);
            colorbar;
            xlabel('Model vocalization');
            ylabel('Actual vocalization');
            title('Confusion Matrix');
            set(gca(), 'Ytick', 1:NstimType);
            set(gca(), 'YTickLabel', StimType);
            set(gca(), 'Xtick', 1:NstimType);
            set(gca(), 'XTickLabel', StimType);
        
             % calculating p(y/x) pour les articles!
            confusion_matrix_CT_p1=zeros(size(confusion_matrix_CallType));
            px=sum(confusion_matrix_CallType, 2);
            for cc=1:size(confusion_matrix_CallType,2)
                confusion_matrix_CT_p1(:,cc)=confusion_matrix_CallType(:,cc)./px;
            end
        
            figure(6);
            imagesc(confusion_matrix_CT_p1);
            colorbar;
            xlabel('Model vocalization');
            ylabel('Actual vocalization');
            title('Confusion Matrix');
            set(gca(), 'Ytick', 1:NstimType);
            set(gca(), 'YTickLabel', StimType);
            set(gca(), 'Xtick', 1:NstimType);
            set(gca(), 'XTickLabel', StimType);
        end
        
        percorrectB(is) = pc;
        mi_confusionB(is) = mi_conf;
        percorrectCT(is) = pcCT;
        mi_confusionCT(is) = mi_confCT;
        
        zdistTB(is) = zdT;
        pzdistTB(is) = pzdT;
        mizdistTB(is) = mi_zdT;
        mizdistncTB(is) = mi_zdT_nc;
        zdistSB(is) = zdS;
        pzdistSB(is) = pzdS;
        mizdistSB(is) = mi_zdS;
        mizdistncSB(is) = mi_zdS_nc;
        confusionMatrix{is} = confusion_matrix;
    end
    

    % Gamma model information calculation
    %onset_time = 200;  % the first 200 ms of the response are removed from the data set
%    fprintf(1, 'Calculating Gamma model information\n');
%    onset_time = 0;
%    spiketrain = spike_times_to_train(NstimType, SpikeTrains, Sections_len, onset_time);
    %[info, noiseentropy, totalentropy, gamma_const, rate_bandwidth, rate_gamma, fano_factor]= gamma_info(spiketrain);

    %gamma_mutual_info = info(1)*1000;                    % Info rates in bits/second
    %gamma_noise_entropy = noiseentropy(1)*1000;       % Noiseentropy in bits/second
%    gamma_spike_rate = mean(mean(spiketrain))*1000;  % Spiking rate is spikes/second

    % Rate Information
    %fprintf(1, 'Calculating Rate Information\n');
    %[rate_info_biased, rate_info_bcorr, rate_info_stderr] = findInfoSR(spiketrain,winSize);

end

%% Fit the neural responses to vocalization type and/or the spectro
fprintf(1, 'fitting the neural responses to vocalization type and/or the spectro\n');
[R2A, ModelPredict, LL, NEC, PValLRatio, h, NeuralResponse, voc, Best_nbPC, Pvalue]=Spectro_Neuro_model(MatfilePath);


%% Store Values
fprintf(1, 'Storing values\n');
Cal.subject=Res.subject;
Cal.StimType=StimType;
Cal.Dsel=Dsel;
Cal.AvgDprime=AvgDprime;
Cal.StdDprime=StdDprime;
Cal.Mean_DP_be=Mean_DP_be;
Cal.Std_DP_be=Std_DP_be;
Cal.Mean_DP_wi=Mean_DP_wi;
Cal.Std_DP_wi=Std_DP_wi;
Cal.Comptype=Comptype;
Cal.AvgRcallcat=AvgRcallcat;
Cal.StdRcallcat=StdRcallcat;
Cal.SSIcat=SSIcat;
Cal.LRIcat=LRIcat;
Dprime=struct();
Dprime.values=dprime;
Dprime.Comptype=AllComptype;
Cal.dprime=Dprime;
Cal.winSize=winSize;
Cal.percorrectB = percorrectB;
Cal.mi_confusionB = mi_confusionB;
Cal.percorrectCT = percorrectCT;
Cal.mi_confusionCT = mi_confusionCT;
Cal.zdistSB = zdistSB;              % Recenters "other" distrubtions to zero to minimize additive variance across songs
Cal.pzdistSB = pzdistSB;             % and collapses distances across all trials of a single stimulus
Cal.mizdistSB = mizdistSB;
Cal.mizdistncSB = mizdistncSB;
%Cal.gamma_mutual_info = gamma_mutual_info;
%Cal.gamma_noise_entropy = gamma_noise_entropy;
Cal.gamma_spike_rate = gamma_spike_rate;
%Cal.totalentropy=totalentropy;
%Cal.gamma_const=gamma_const;
%Cal.rate_bandwidth=rate_bandwidth;
%Cal.rate_gamma=rate_gamma;
%Cal.fano_factor=fano_factor;
%Cal.rate_info_biased=rate_info_biased;
%Cal.rate_info_bcorr=rate_info_bcorr;
%Cal.rate_info_stderr=rate_info_stderr;
Cal.Best_nbPC=Best_nbPC;
Cal.R2A=R2A;
Cal.ModelPredict = ModelPredict;
Cal.LogLikelihood = LL;
Cal.Pvalue=Pvalue; %result of the anova on each model
Cal.PValLRatio = PValLRatio;
Cal.SignifModelCompare = h;
Cal.NeuralResponse = NeuralResponse;
Cal.VocType=voc;
%Cal.confusionMatrix = confusionMatrix;

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
clear Res Dprime Cal  dsel AvgDprime StdDprime Comptype AvgRcallcat StdRcallcat winSize percorrectB mi_confusionB zdistSB pzdistSB mizdistSB mizdistncSB gamma_spike_rate R2A confusionMatrix
%clear gamma_mutual_info gamma_noise_entropy totalentropy gamma_const
%rate_bandwidth rate_gamma fano_factor rate_info_biased rate_info_bcorr rate_info_stderr
end 
