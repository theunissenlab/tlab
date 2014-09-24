function [R2A, ModelPredict, LL, NEC, PValLRatio, h, y, voc, Best_nbPC, Pvalue]=Spectro_Neuro_model(MatfilePath)
FIG = 0;
%nPC=100;
NBPC=1:10:110; % to see how adjusted R2 evolve with NB of PC
%NBPC=38; %Mean Nb of PC that gives the best results of adjusted R2 as observed by running this code with the previous line with the script optimalPC_ModelNeuralResponse.m on 10/02/2012
R2A=zeros(3,1);
ModelPredict = cell(3,1);
LL= zeros(3,1);
Pvalue=zeros(3,1);              %results from the anova on each model
NEC = zeros(3,1);
PValLRatio = zeros(3,1);
h = zeros(3,1);
R2A_temp=zeros(3,length(NBPC));
ModelPredict_temp = cell(3,length(NBPC));
LL_temp= zeros(3,length(NBPC));
Pvalue_temp=zeros(3,length(NBPC));
NEC_temp = zeros(3,length(NBPC));

%% Load the unit matfile
%Res=load('/Users/elie/Documents/MATLAB/data/matfile/GreBlu9508M/ZS_Site2_L1100R1450_21.mat')
Res = load(MatfilePath);
DataSel=zeros(1,length(Res.VocType));
nvoc=1;
voctype=Res.VocType;
for dd=1:length(voctype);
    if strcmp(voctype{dd}, 'Ag')
        DataSel(nvoc)=dd;
        nvoc=nvoc+1;
    elseif strcmp(voctype{dd}, 'Be')
        DataSel(nvoc)=dd;
        nvoc=nvoc+1;
    elseif strcmp(voctype{dd}, 'DC')
        DataSel(nvoc)=dd;
        nvoc=nvoc+1;
    elseif strcmp(voctype{dd}, 'Di')
        DataSel(nvoc)=dd;
        nvoc=nvoc+1;
    elseif strcmp(voctype{dd}, 'LT')
        DataSel(nvoc)=dd;
        nvoc=nvoc+1;
    elseif strcmp(voctype{dd}, 'Ne')
        DataSel(nvoc)=dd;
        nvoc=nvoc+1;
    elseif strcmp(voctype{dd}, 'Te')
        DataSel(nvoc)=dd;
        nvoc=nvoc+1;
    elseif strcmp(voctype{dd}, 'Th')
        DataSel(nvoc)=dd;
        nvoc=nvoc+1;
    elseif strcmp(voctype{dd}, 'song')
        DataSel(nvoc)=dd;
        nvoc=nvoc+1;
    elseif strcmp(voctype{dd}, 'Wh')
        DataSel(nvoc)=dd;
        nvoc=nvoc+1;
    end
end
if (nvoc-1)==length(Res.MeanRate)
    y= cell2mat(Res.MeanRate);
    x = cell2mat(Res.Spectro);
    voc = Res.VocType;
else
    DataSel=DataSel(1:(nvoc-1));
    y= cell2mat(Res.MeanRate(DataSel));
    x = cell2mat(Res.Spectro(DataSel));
    voc = Res.VocType(DataSel);
end
x = 20*log10(abs(x)); % here we take the log of the spectro but pb: x
MAXI = max(max(x));
x(find(x<(MAXI-80)))=MAXI-80;
%contains 0 and log(x) gives inf and princomp does not support inf... :-(

%% Calculate the principal components of the spectro
fprintf(1, 'Calculate PC of spectro\n');
[COEFF,SCORE,latent,tsquare]=princomp(x,'econ');
if FIG==1
    figure(1)
    plot(cumsum(latent/sum(latent)))
end
pValue=0;
jj=0;
%while jj<length(NBPC)
while pValue<0.05 && jj<length(NBPC)
    jj = jj + 1;
    nPC=NBPC(jj);
    %% Fit the neural response of the unit with the first 100 PC of the spectro
    % and calculate the STRF
    %fprintf(1, 'Calculating STRF using the %d first PC of the spectro\n', nPC);
    %mdl=LinearModel.fit(SCORE(:,1:nPC),y);
    %PCSTRF=mdl.Coefficients.Estimate(2:end);
    %STRF=COEFF(:,1:nPC)*PCSTRF;
    %STRFM=reshape(STRF,length(Res.Spectrofo), length(Res.Spectroto));
    if FIG==1
        figure(jj+1)
        imagesc(Res.Spectroto, Res.Spectrofo, STRFM)
        axis xy
        pause
    end

    %% Fit the neural data with the PC of the spectro and/or the vocalization type and retrieve RMSE 
    %ds=dataset(SCORE(:,1:100),y);
    fprintf('Constructing three models of the neural response (spike rate) with the following variables:\n PC of spectro\nVocalization type\nPC of spectro + vocalization type\n');
    ds=dataset();

    for ii=1:nPC
        ds.(sprintf('SCORE%d',ii)) = SCORE(:,ii);
    end

    ds2=dataset();
    ds2.Voctype=ordinal(voc);
    ds2.y=y;
    
    % Plot of the predicted spike rate given the voctype
    if FIG==1
        CoeffEstimates=mdl2.Coefficients.Estimate(1:end);
        MeanValues=CoeffEstimates + [0 ; repmat(CoeffEstimates(1), (length(CoeffEstimates)-1),1)];
        plot(MeanValues, 1:length(MeanValues));
        pause
    end
    
    %keep on going with other models
    ds3=ds;
    ds.y=y;
    ds3.Voctype=ordinal(voc);
    ds3.y=y;

    mdl=LinearModel.fit(ds);    %Model with PC of spectro only
    R2A_temp(1,jj)=mdl.Rsquared.Adjusted;
    ModelPredict_temp{1,jj}=mdl.predict;
    LL_temp(1,jj) = mdl.LogLikelihood;
    NEC_temp(1,jj) = mdl.NumEstimatedCoefficients;
    tbl=anova(mdl,'summary');
    Pvalue_temp(1,jj)=tbl.pValue(2);
    mdl2=LinearModel.fit(ds2);  %Model with  VocType only
    R2A_temp(2,jj)=mdl2.Rsquared.Adjusted;
    ModelPredict_temp{2,jj}=mdl2.predict;
    LL_temp(2,jj) = mdl2.LogLikelihood;
    NEC_temp(2,jj) = mdl2.NumEstimatedCoefficients;
    tbl2=anova(mdl2,'summary');
    Pvalue_temp(2,jj)=tbl2.pValue(2);
    mdl3=LinearModel.fit(ds3);  %Model with both PC of spectro and VocType
    R2A_temp(3,jj)=mdl3.Rsquared.Adjusted;
    ModelPredict_temp{3,jj}=mdl3.predict;
    LL_temp(3,jj) = mdl3.LogLikelihood;
    NEC_temp(3,jj) = mdl3.NumEstimatedCoefficients;
    tbl3=anova(mdl3,'summary');
    Pvalue_temp(3,jj)=tbl3.pValue(2);
    if jj>1
        dof = NEC_temp(1,jj) - NEC_temp(1,(jj-1));
        [h_temp,pValue,stat,cValue] = lratiotest(LL_temp(1,(jj)),LL_temp(1,jj-1),dof);
    end
    pValue
end
Best_nbPC=NBPC(jj-1); %Best number of PC calculated on the acoustic model and use for this unit for the calculations of the 3 models
fprintf('best nb of PC= %d\n', Best_nbPC);
R2A(1:3,1) = R2A_temp(1:3,(jj-1));
ModelPredict{1,1} = ModelPredict_temp{1,(jj-1)};
ModelPredict{2,1} = ModelPredict_temp{2,(jj-1)};
ModelPredict{3,1} = ModelPredict_temp{3,(jj-1)};
LL(1:3,1) = LL_temp(1:3,(jj-1));
NEC(1:3,1) =NEC_temp(1:3,(jj-1));
Pvalue(1:3,1)=Pvalue_temp(1:3,(jj-1));
if NEC(1,1)>NEC(2,1)
    [h(1),PValLRatio(1,1),stat,cValue] = lratiotest(LL(1,1),LL(2,1),NEC(1,1) - NEC(2,1));
elseif NEC(1,1)<NEC(2,1)
    [h(1),PValLRatio(1,1),stat,cValue] = lratiotest(LL(2,1),LL(1,1),NEC(2,1) - NEC(1,1));
else
    fprintf('same degree of freedom for Lratiotest between Acoustic and semantic');
    [h(1),PValLRatio(1,1),stat,cValue] = lratiotest(LL(1,1),LL(2,1),NEC(1,1) - NEC(2,1));
end
[h(2),PValLRatio(2,1),stat,cValue] = lratiotest(LL(3,1),LL(1,1),NEC(3,1) - NEC(1,1));
[h(3),PValLRatio(3,1),stat,cValue] = lratiotest(LL(3,1),LL(2,1),NEC(3,1) - NEC(2,1));

%fprintf('the RMSE of the models are:\n%d first PC of spectro: %f\nVocType only: %f\n%d first PC of spectro + VocType: %f\n', nPC, RMSE(1), RMSE(2), nPC, RMSE(3));
if FIG==2
    figure(12)
    plot(NBPC, R2A_temp(1,:), 'rs-', NBPC, R2A_temp(2,:), 'co-', NBPC, R2A_temp(3,:), 'g*-');
    axis([0 100 0 1])
end
end

