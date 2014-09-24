%% Perform a regularized DFA based on Spectrogrammes for all the zebra finch vocalizations
% This is the version of DFA_Calls_Julie that does the individual
% recognition.

% Read the data base produced by VocSectioningAmp.m 
load /Users/frederictheunissen/Documents/Data/Julie/analysis/vocCuts.mat

% To make RAM space, clear the sound files
clear soundCutsTot;
pack;

% Set a 100 dB range threshold
maxAmp = max(max(spectroCutsTot));
minAmp = maxAmp - 100;
spectroCutsTot(spectroCutsTot<minAmp) = minAmp;


%% Perform a regularized PCA
% At this point this is just the vanilla PCA

[Coeff, Score, latent] = princomp(spectroCutsTot, 'econ'); 
% Here Coeff are the PCs organized in columns and Score is the projection
% of the spectrogram on the PCs. One could get the Scores from the Coeffs.
% mSpectro = mean(spectroCutsTot);
% xSpectro = spectroCutsTot - repmat(mSpectro, size(spectroCutsTot,1), 1);
% Score = Spectro*Coeff

% We only save nb Scores and clear the spectrograms to make space
clear spectroCutsTot;
pack;
save /Users/frederictheunissen/Documents/Data/Julie/analysis/vocPCA.mat

%% Start here if you already done the PCA
load /Users/frederictheunissen/Documents/Data/Julie/analysis/vocPCA.mat


%% Some quick clean up
% Remove unknown from data set
unknownInd = find(strncmp(birdNameCuts,'Unknown', 7));
ScoreInd = Score;
ScoreInd(unknownInd,:) = [];
birdNameCutsInd = birdNameCuts;
birdNameCutsInd(unknownInd) = [];
vocTypeCutsInd = vocTypeCuts;
vocTypeCutsInd(unknownInd) = [];
ncutsTotInd = length(birdNameCutsInd);


%%  Perform the classical DFA on all of the Data to obtain DFA for Identity
% Here we look at the effect of the number of PCs on the DFAs

% Calculate first for 300 PCs
nb = 300;
[nDF, p, stats] = manova1(ScoreInd(:, 1:nb), birdNameCutsInd);
[mean_grp, sem_grp, meanCI_grp, range_grp, name_grp] = grpstats(stats.canon(:,1:nDF),birdNameCutsInd', {'mean', 'sem', 'meanci', 'range', 'gname'});
ngroups = size(mean_grp,1);
PC_DF_300 = Coeff(:, 1:nb) * stats.eigenvec(:, 1:nDF);

DF_corr = zeros(5, nDF);

for nbi=1:5
    nb = nbi*50;   % number of PCs
    [nDF, p, stats] = manova1(ScoreInd(:, 1:nb), birdNameCutsInd);
    [mean_grp, sem_grp, meanCI_grp, range_grp, name_grp] = grpstats(stats.canon(:,1:nDF),birdNameCutsInd', {'mean', 'sem', 'meanci', 'range', 'gname'});
    ngroups = size(mean_grp,1);
    
    %  Display the significant DFA in the spectrogram space
    PC_DF = Coeff(:, 1:nb) * stats.eigenvec(:, 1:nDF);
    
    for i = 1:nDF
        DF_corr(nbi,i) = corr(PC_DF(:,i), PC_DF_300(:,i));
    end
    
    % Find color scale
    clear cmin cmax clims
    cmin = min(min(PC_DF(:, 1:nDF)));
    cmax = max(max(PC_DF(:, 1:nDF)));
    cabs = max(abs(cmin), abs(cmax));
    clims = [-cabs cabs];
    nf = length(fo);
    nt = length(to);
    
    % Plot significant DFAs
    figure(3);
    for i=1:nDF
        subplot(6,nDF,i+(nbi-1)*nDF);
        PC_Spect = reshape(PC_DF(:,i), nf, nt);
        imagesc(to, fo, PC_Spect, clims)
        if i == 1
            title(sprintf('Indiv DFs using %d PC)', nb));
        end
        axis xy;
        axis([0 0.2 0 8000]);
    end
end

% Plot the last line
for i=1:nDF
        subplot(6,nDF,i+5*nDF);
        PC_Spect = reshape(PC_DF_300(:,i), nf, nt);
        imagesc(to, fo, PC_Spect, clims)
        if i == 1
            title(sprintf('Indv DFs using %d PC)', 300));
        end
        axis xy;
        axis([0 0.2 0 8000]);
 end

% Plot the coefficients of the DFA for n=50, n=100, n=200 PCs
figure(4);
% Using 50 PCs
subplot(1,4,1);
[nDF, p, stats] = manova1(ScoreInd(:, 1:50), birdNameCutsInd); 
weigthsPC_DFA = zeros(300, nDF);
weigthsPC_DFA(1:50, :) = abs(stats.eigenvec(1:50, 1:nDF));
imagesc(1:nDF, 1:300, weigthsPC_DFA);
ylabel('PCA Vector');
title('50 PCs');

% Using 100 PCs
subplot(1,4,2);
[nDF, p, stats] = manova1(ScoreInd(:, 1:100), birdNameCutsInd); 
weigthsPC_DFA = zeros(300, nDF);
weigthsPC_DFA(1:100, :) = abs(stats.eigenvec(1:100, 1:nDF));
imagesc(1:nDF, 1:300, weigthsPC_DFA);
ylabel('PCA Vector');
title('100 PCs');

% Using 200 PCs
subplot(1,4,3);
[nDF, p, stats] = manova1(ScoreInd(:, 1:200), birdNameCutsInd); 
weigthsPC_DFA = zeros(300, nDF);
weigthsPC_DFA(1:200, :) = abs(stats.eigenvec(1:200, 1:nDF));
imagesc(1:nDF, 1:300, weigthsPC_DFA);
ylabel('PCA Vector');
title('200 PCs');

% Using 300 PCs
subplot(1,4,4);
[nDF, p, stats] = manova1(ScoreInd(:, 1:300), birdNameCutsInd); 
weigthsPC_DFA = zeros(300, nDF);
weigthsPC_DFA(1:300, :) = abs(stats.eigenvec(1:300, 1:nDF));
imagesc(1:nDF, 1:300, weigthsPC_DFA);
ylabel('PCA Vector');
title('300 PCs');

% Plot the correlation coefficient between DFA's and the one calculated
% with 300 PCs.
figure(5);
plot(50:50:250, abs(DF_corr));
xlabel('Number of PCs');
ylabel('Correlation with DFA at 300');
legend('DFA1', 'DFA2', 'DFA3', 'DFA4', 'DFA5', 'DFA6', 'DFA7', 'DFA8');


%%  Perform a boot-strap to calculate confusion matrix and percent correct
% The boot-strap is also performed for various DFA to obtain the optimal
% number.

n_boot = 1000;     % Number of permutations
n_valid = fix(0.1*ncutsTotInd); % Save 10% of the data for cross-validation

PCC_info = struct('nb', {50, 100, 150, 200, 250, 300}, 'PCC_Total', 0, 'PCC_M', 0, 'PCC_group', zeros(ngroups, ngroups));

for inb = 1:length(PCC_info)
    
    % Number of PCA used in the DFA
    nb = PCC_info(inb).nb;
    
    % Allocate space for distance vector and confusion matrix
    Dist = zeros(1, ngroups);
    ConfMat = zeros(ngroups, ngroups);
    
    for iboot=1:n_boot
        
        ind_valid = randsample(ncutsTotInd, n_valid);    % index of the validation calls
        
        % Separate data into fitting and validation
        X_valid = ScoreInd(ind_valid, 1:nb);
        X_fit = ScoreInd(:, 1:nb);
        X_fit(ind_valid, :) = [];
        
        % Similarly for the group labels.
        Group_valid = birdNameCutsInd(ind_valid);
        Group_fit = birdNameCutsInd;
        Group_fit(ind_valid) = [];
        
        % Perform the linear DFA using manova1 for the training set
        [nDF, p, stats] = manova1(X_fit, Group_fit);
        [mean_bgrp, sem_bgrp, meanbCI_grp, range_bgrp, name_bgrp] = grpstats(stats.canon(:,1:nDF),Group_fit', {'mean', 'sem', 'meanci', 'range', 'gname'});
        nbgroups = size(mean_bgrp,1);
        
        % Project the validation data set into the DFA.
        mean_X_fit = mean(X_fit);
        Xc = X_valid - repmat(mean_X_fit, size(X_valid,1), 1);
        Canon = Xc*stats.eigenvec(:, 1:nDF);
        
        % Use Euclidian Distances
        for i = 1:n_valid
            for j = 1:nbgroups
                Dist(j) = sqrt((Canon(i,:) - mean_bgrp(j,:))*(Canon(i,:) - mean_bgrp(j,:))');
                if strcmp(name_bgrp(j),Group_valid(i))
                    k_actual = j;
                end
            end
            k_guess = find(Dist == min(Dist), 1, 'first');
            
            % Just in case a group is missing find the index that corresponds
            % to the groups when all the data is taken into account.
            for j=1:ngroups
                if strcmp(name_grp(j), name_bgrp(k_actual))
                    k_actual_all = j;
                    break;
                end
            end
            for j=1:ngroups
                if strcmp(name_grp(j), name_bgrp(k_guess))
                    k_guess_all = j;
                    break;
                end
            end
            
            ConfMat(k_actual_all, k_guess_all) = ConfMat(k_actual_all, k_guess_all) + 1;
        end
        
        
        
    end
    
    PCC_Total = 100.0*sum(diag(ConfMat))./(n_valid*n_boot);
    PCC_group = zeros(ngroups, ngroups);
    for i = 1:ngroups
        for j = 1:ngroups
            PCC_group(i,j) = ConfMat(i,j) / sum(ConfMat(i, :), 2) * 100; % sum(.., 2) = somme par ligne
        end
    end
    PCC_M = mean(diag(PCC_group));
    
    % Store the information
    
    PCC_info(inb).PCC_Total = PCC_Total;
    PCC_info(inb).PCC_M = PCC_M;
    PCC_info(inb).PCC_group = PCC_group;
    
end
% Display confusion Matrix

figure(6);
for inb = 1:length(PCC_info)
    subplot(1,length(PCC_info),inb);
    imagesc(PCC_info(inb).PCC_group);
    xlabel('Guess');
    ylabel('Actual');
    colormap(gray);
    colorbar;
    title(sprintf('Confusion Matrix RF %.1f%%(%.1f%%) Correct', PCC_info(inb).PCC_Total, PCC_info(inb).PCC_M));
    set(gca(), 'Ytick', 1:ngroups);
    set(gca(), 'YTickLabel', name_grp);
    set(gca(), 'Xtick', 1:ngroups);
    set(gca(), 'XTickLabel', name_grp);
end

figure(7);
plot([PCC_info.nb], [PCC_info.PCC_Total],'k');
hold on;
plot([PCC_info.nb], [PCC_info.PCC_M],'r');
xlabel('Number of PCs');
ylabel('Correct Classification');
legend('Total', 'Averaged per Ind');

% save the DFA
save vocIndDFA.mat fo to PC_DF PCC_info ngroups name_grp birdNameCutsInd ncutsTotInd ScoreInd vocTypeCutsInd

%% From these results we choose nb=150 PCs and make nice plots + save some data
% You can load vocTypeDFA.mat and vocPCA.mat and start here to make these
% final plots
% load /Users/frederictheunissen/Documents/Data/Julie/analysis/vocPCA.mat
% load /Users/frederictheunissen/Documents/Data/Julie/analysis/vocIndDFA.mat

% Find the birds that are juveniles
% These are the indices of all the Ne and LT calls
ind = find(strcmp([vocTypeCutsInd],'Be')|strcmp([vocTypeCutsInd],'LT'));
juvieNames = unique([birdNameCutsInd(ind)]);

for i=1:length(name_grp)
    juvieFlg = 0;
    for j=1:length(juvieNames)
        if strcmp(name_grp(i), juvieNames(j))
            juvieFlg = 1;
            break;
        end
    end
    
    if juvieFlg
        name_grpA{i} = [name_grp{i},'J'];
    else
        name_grpA{i} = [name_grp{i},'A'];
    end
end


% A nice confusion matrix.
figure(8);
inb = 3;   % this corresponds to 150 pcs

imagesc(PCC_info(inb).PCC_group);
xlabel('Guess');
ylabel('Actual');
colormap(gray);
colorbar;
title(sprintf('nPCA = %d PC Total= %.1f%% PCC Mean = %.1f%%', PCC_info(inb).nb, PCC_info(inb).PCC_Total, PCC_info(inb).PCC_M));
set(gca(), 'Ytick', 1:ngroups);
set(gca(), 'YTickLabel', name_grpA);
set(gca(), 'Xtick', 1:ngroups);
set(gca(), 'XTickLabel', name_grpA);

figure(9);
nb = 150;
[nDF, p, stats] = manova1(ScoreInd(:, 1:nb), birdNameCutsInd);
[mean_grp, std_grp, name_grp] = grpstats(stats.canon(:,1:nDF),birdNameCutsInd', {'mean', 'std', 'gname'});

    
%  Display the significant DFA in the spectrogram space
PC_DF = Coeff(:, 1:nb) * stats.eigenvec(:, 1:nDF);

% Find color scale
clear cmin cmax clims
cmin = min(min(PC_DF(:, 1:nDF)));
cmax = max(max(PC_DF(:, 1:nDF)));
cabs = max(abs(cmin), abs(cmax));
clims = [-cabs cabs];
nf = length(fo);
nt = length(to);


for i=1:nDF
    subplot(1,nDF,i);
    PC_Spect = reshape(PC_DF(:,i), nf, nt);
    imagesc(to, fo, PC_Spect, clims)
    
    title(sprintf('DFA%d', i));
    
    axis xy;
    axis([0 0.2 0 8000]);
    
    if (i == nDF)
        xlabel('Time (s)');
    end
    if (i == 1)
        ylabel('Frequency (Hz)');
    end
end

% Make DF1 vs DFX plots
figure(10);
colorVals = colormap('lines');
ncolor = length(colorVals);

for iDF=2:nDF
    subplot(2,ceil((nDF-1)/2), iDF-1);
    for ig=1:ngroups
        
        plot(mean_grp(ig,1), mean_grp(ig,iDF), 's', 'MarkerSize', 10, ...
                'color', colorVals(mod(ig,ncolor),:),'MarkerEdgeColor','k',...
                'MarkerFaceColor',colorVals(mod(ig,ncolor),:));
        hold on;
    end
    xlabel('DF1');
    ylabel(sprintf('DF%d', iDF));
    axis([-8 8 -8 8]);
    axis square;
    if iDF == nDF
        legend(name_grpA);
    end
    hold off;
end
    

%% Repeat the classification with the Random Forest Algorithm   

% The random tree does its own bootstrat and we don't need a training and
% validation set.  
nb = 150;                    % Here we choose 150 PCs


Group_names = unique(birdNameCutsInd);
ngroups = length(Group_names);

for i=1:ngroups
    juvieFlg = 0;
    for j=1:length(juvieNames)
        if strcmp(name_grp(i), juvieNames(j))
            juvieFlg = 1;
            break;
        end
    end
    
    if juvieFlg
        Group_namesA{i} = [Group_names{i},'J'];
    else
        Group_namesA{i} = [Group_names{i},'A'];
    end
end



ConfMat = zeros(ngroups, ngroups);   

B = TreeBagger(300, ScoreInd(:, 1:nb), birdNameCutsInd, 'OOBPred', 'on', 'priorprob', 'equal', 'MinLeaf', 5, 'NPrint', 10);

Group_predict = oobPredict(B);   % This returns the predictions for the out of bag values.

n_valid = length(birdNameCutsInd);   % this is the total number of observations

for i = 1:n_valid
    k_actual = find(strcmp(Group_names,birdNameCutsInd(i)));
    k_guess = find(strcmp(Group_names,Group_predict(i)));
    ConfMat(k_actual, k_guess) = ConfMat(k_actual, k_guess) + 1;
end

PCC_Total = 100.0*sum(diag(ConfMat))./(n_valid);

PCC_group = zeros(ngroups, ngroups);
for i = 1:ngroups
    for j = 1:ngroups
        PCC_group(i,j) = ConfMat(i,j) / sum(ConfMat(i, :), 2) * 100; 
    end
end
PCC_M = mean(diag(PCC_group));

% Plot out of bag error
figure(11);
plot(1-oobError(B));
xlabel('number of grown trees')
ylabel('out-of-bag Prediction')
title('Random Forest Performance for Call Classification');

% Display confusion Matrix
figure(12);
imagesc(PCC_group);
xlabel('Guess');
ylabel('Actual');
colormap(gray);
colorbar;
title(sprintf('Confusion Matrix RF %.1f%%(%.1f%%) Correct', PCC_Total, PCC_M));
set(gca(), 'Ytick', 1:ngroups);
set(gca(), 'YTickLabel', Group_namesA);
set(gca(), 'Xtick', 1:ngroups);
set(gca(), 'XTickLabel', Group_namesA);

%% Repeat the calcultion but separating juveniles from adults and keeping track of vocalization type.

indJuvie = find(strcmp([vocTypeCutsInd],'Be')|strcmp([vocTypeCutsInd],'LT'));
juvieNames = unique([birdNameCutsInd(indJuvie)]);

indAdult = find(~strcmp([vocTypeCutsInd],'Be')&~strcmp([vocTypeCutsInd],'LT'));
adultNames = unique([birdNameCutsInd(indAdult)]);

%% First Juveniles


n_boot = 1000;     % Number of permutations


ScoreJuv = ScoreInd(indJuvie,:);
vocTypeJuv = vocTypeCutsInd(indJuvie);
birdNameJuv = birdNameCutsInd(indJuvie);
juvieVocs = unique(vocTypeJuv);
ncallsJuv = length(indJuvie);


n_valid = fix(0.1*ncallsJuv); % Save 10% of the data for cross-validation


% Number of PCA used in the DFA
nb = 150;

% Allocate space for distance vector and confusion matrix
nInd = length(juvieNames);
nVoc = length(juvieVocs);

Dist = zeros(1, nInd);
ConfMat = zeros(nVoc, nInd, nInd);

for iboot=1:n_boot
    
    ind_valid = randsample(ncallsJuv, n_valid);    % index of the validation calls
    
    % Separate data into fitting and validation
    X_valid = ScoreJuv(ind_valid, 1:nb);
    X_fit = ScoreJuv(:, 1:nb);
    X_fit(ind_valid, :) = [];
    
    % Similarly for the group labels.
    Group_valid = birdNameJuv(ind_valid);
    Group_fit = birdNameJuv;
    Group_fit(ind_valid) = [];
    
    % And for voc type labels
    Type_valid = vocTypeJuv(ind_valid);
    Type_fit = vocTypeJuv;
    Type_fit(ind_valid) = [];
    
    
    % Perform the linear DFA using manova1 for the training set
    [nDF, p, stats] = manova1(X_fit, Group_fit);
    [mean_bgrp, sem_bgrp, meanbCI_grp, range_bgrp, name_bgrp] = grpstats(stats.canon(:,1:nDF),Group_fit', {'mean', 'sem', 'meanci', 'range', 'gname'});
    nbgroups = size(mean_bgrp,1);
    
    % Project the validation data set into the DFA.
    mean_X_fit = mean(X_fit);
    Xc = X_valid - repmat(mean_X_fit, size(X_valid,1), 1);
    Canon = Xc*stats.eigenvec(:, 1:nDF);
    
    % Use Euclidian Distances
    for i = 1:n_valid
        iVocs = find(strcmp([juvieVocs], Type_valid(i))); 
        for j = 1:nbgroups
            Dist(j) = sqrt((Canon(i,:) - mean_bgrp(j,:))*(Canon(i,:) - mean_bgrp(j,:))');
            if strcmp(name_bgrp(j),Group_valid(i))
                k_actual = j;
            end
        end
        k_guess = find(Dist == min(Dist), 1, 'first');
        
        % Just in case a group is missing find the index that corresponds
        % to the groups when all the data is taken into account.
        for j=1:nInd
            if strcmp(juvieNames(j), name_bgrp(k_actual))
                k_actual_all = j;
                break;
            end
        end
        for j=1:nInd
            if strcmp(juvieNames(j), name_bgrp(k_guess))
                k_guess_all = j;
                break;
            end
        end
        
        ConfMat(iVocs, k_actual_all, k_guess_all) = ConfMat(iVocs, k_actual_all, k_guess_all) + 1;
    end
    
        
end

% First show the overall Performance
ConfMatAll = squeeze(sum(ConfMat, 1));
PCC_Total = 100.0*sum(diag(ConfMatAll))./(n_valid*n_boot);
PCC_group = zeros(nInd, nInd);
for i = 1:nInd
    for j = 1:nInd
        PCC_group(i,j) = ConfMatAll(i,j) / sum(ConfMatAll(i, :), 2) * 100; % sum(.., 2) = somme par ligne
    end
end
PCC_M = mean(diag(PCC_group));

% A nice confusion matrix.
figure(13);

imagesc(PCC_group);
xlabel('Guess');
ylabel('Actual');
colormap(gray);
colorbar;
title(sprintf('Total= %.1f%% PCC Mean = %.1f%%',PCC_Total, PCC_M));
set(gca(), 'Ytick', 1:nInd);
set(gca(), 'YTickLabel', juvieNames);
set(gca(), 'Xtick', 1:nInd);
set(gca(), 'XTickLabel', juvieNames);

% Now per vocalization type
figure(14);
for iVoc=1:nVoc
    subplot(1,nVoc, iVoc);
    ConfMatVoc = squeeze(ConfMat(iVoc,:,:));
    n_tot = sum(sum(ConfMatVoc));
    PCC_Total = 100.0*sum(diag(ConfMatVoc))./n_tot;
    PCC_group = zeros(nInd, nInd);
    for i = 1:nInd
        for j = 1:nInd
            PCC_group(i,j) = ConfMatVoc(i,j) / sum(ConfMatVoc(i, :), 2) * 100; % sum(.., 2) = somme par ligne
        end
    end
    PCC_M = mean(diag(PCC_group));

    imagesc(PCC_group);
    xlabel('Guess');
    ylabel('Actual');
    colormap(gray);
    caxis([0 100]);
    if iVoc == nVoc
        colorbar;
    end
    
    % Find the number for each individual
    numberInd = zeros(1, nInd);
    for iInd=1:nInd
        numberInd(iInd) = length(find(strcmp([birdNameJuv], juvieNames(iInd)) & strcmp([vocTypeJuv], juvieVocs(iVoc))));
        juvieNamesNum{iInd} = sprintf('%s (%d)', juvieNames{iInd}, numberInd(iInd));
    end
    
    title(sprintf('%s Total= %.1f%% PCC Mean = %.1f%%', juvieVocs{iVoc}, PCC_Total, PCC_M));
    set(gca(), 'Ytick', 1:nInd);
    set(gca(), 'YTickLabel', juvieNamesNum);
    set(gca(), 'Xtick', 1:nInd);
    set(gca(), 'XTickLabel', juvieNames);
end
    
%% Now Adults


n_boot = 1000;     % Number of permutations


ScoreAd = ScoreInd(indAdult,:);
vocTypeAd = vocTypeCutsInd(indAdult);
birdNameAd = birdNameCutsInd(indAdult);
adultVocs = unique(vocTypeAd);
ncallsAd = length(indAdult);


n_valid = fix(0.1*ncallsAd); % Save 10% of the data for cross-validation


% Number of PCA used in the DFA
nb = 150;

% Allocate space for distance vector and confusion matrix
nInd = length(adultNames);
nVoc = length(adultVocs);

Dist = zeros(1, nInd);
ConfMat = zeros(nVoc, nInd, nInd);

for iboot=1:n_boot
    
    ind_valid = randsample(ncallsAd, n_valid);    % index of the validation calls
    
    % Separate data into fitting and validation
    X_valid = ScoreAd(ind_valid, 1:nb);
    X_fit = ScoreAd(:, 1:nb);
    X_fit(ind_valid, :) = [];
    
    % Similarly for the group labels.
    Group_valid = birdNameAd(ind_valid);
    Group_fit = birdNameAd;
    Group_fit(ind_valid) = [];
    
    % And for voc type labels
    Type_valid = vocTypeAd(ind_valid);
    Type_fit = vocTypeAd;
    Type_fit(ind_valid) = [];
    
    
    % Perform the linear DFA using manova1 for the training set
    [nDF, p, stats] = manova1(X_fit, Group_fit);
    [mean_bgrp, sem_bgrp, meanbCI_grp, range_bgrp, name_bgrp] = grpstats(stats.canon(:,1:nDF),Group_fit', {'mean', 'sem', 'meanci', 'range', 'gname'});
    nbgroups = size(mean_bgrp,1);
    
    % Project the validation data set into the DFA.
    mean_X_fit = mean(X_fit);
    Xc = X_valid - repmat(mean_X_fit, size(X_valid,1), 1);
    Canon = Xc*stats.eigenvec(:, 1:nDF);
    
    % Use Euclidian Distances
    for i = 1:n_valid
        iVocs = find(strcmp([adultVocs], Type_valid(i))); 
        for j = 1:nbgroups
            Dist(j) = sqrt((Canon(i,:) - mean_bgrp(j,:))*(Canon(i,:) - mean_bgrp(j,:))');
            if strcmp(name_bgrp(j),Group_valid(i))
                k_actual = j;
            end
        end
        k_guess = find(Dist == min(Dist), 1, 'first');
        
        % Just in case a group is missing find the index that corresponds
        % to the groups when all the data is taken into account.
        for j=1:nInd
            if strcmp(adultNames(j), name_bgrp(k_actual))
                k_actual_all = j;
                break;
            end
        end
        for j=1:nInd
            if strcmp(adultNames(j), name_bgrp(k_guess))
                k_guess_all = j;
                break;
            end
        end
        
        ConfMat(iVocs, k_actual_all, k_guess_all) = ConfMat(iVocs, k_actual_all, k_guess_all) + 1;
    end
    
        
end

% First show the overall Performance
ConfMatAll = squeeze(sum(ConfMat, 1));
PCC_Total = 100.0*sum(diag(ConfMatAll))./(n_valid*n_boot);
PCC_group = zeros(nInd, nInd);
for i = 1:nInd
    for j = 1:nInd
        PCC_group(i,j) = ConfMatAll(i,j) / sum(ConfMatAll(i, :), 2) * 100; % sum(.., 2) = somme par ligne
    end
end
PCC_M = mean(diag(PCC_group));

% A nice confusion matrix.
figure(15);

imagesc(PCC_group);
xlabel('Guess');
ylabel('Actual');
colormap(gray);
colorbar;
title(sprintf('Total= %.1f%% PCC Mean = %.1f%%',PCC_Total, PCC_M));
set(gca(), 'Ytick', 1:nInd);
set(gca(), 'YTickLabel', adultNames);
set(gca(), 'Xtick', 1:nInd);
set(gca(), 'XTickLabel', adultNames);

% Now per vocalization type
figure(16);
for iVoc=1:nVoc
    subplot(3,ceil(nVoc/3), iVoc);
    ConfMatVoc = squeeze(ConfMat(iVoc,:,:));
    n_tot = sum(sum(ConfMatVoc));
    PCC_Total = 100.0*sum(diag(ConfMatVoc))./n_tot;
    PCC_group = zeros(nInd, nInd);
    for i = 1:nInd
        for j = 1:nInd
            PCC_group(i,j) = ConfMatVoc(i,j) / sum(ConfMatVoc(i, :), 2) * 100; % sum(.., 2) = somme par ligne
        end
    end
    PCC_M = mean(diag(PCC_group));

    imagesc(PCC_group);
    xlabel('Guess');
    ylabel('Actual');
    colormap(gray);
    caxis([0 100]);
%     if iVoc == nVoc
%         colorbar;
%     end
    
    % Find the number for each individual
    numberInd = zeros(1, nInd);
    for iInd=1:nInd
        numberInd(iInd) = length(find(strcmp([birdNameAd], adultNames(iInd)) & strcmp([vocTypeAd], adultVocs(iVoc))));
        if iVoc == 1
            adultNamesNum{iInd} = sprintf('%s (%d)', adultNames{iInd}, numberInd(iInd));
        else
            adultNamesNum{iInd} = sprintf('(%d)', numberInd(iInd));
        end
    end
    
    title(sprintf('%s Total= %.1f%% PCC Mean = %.1f%%', adultVocs{iVoc}, PCC_Total, PCC_M));
    set(gca(), 'Ytick', 1:nInd);
    set(gca(), 'YTickLabel', adultNamesNum);
    set(gca(), 'Xtick', 1:nInd);
    % This is commented out because it is too busy with it...
    % set(gca(), 'XTickLabel', adultNames);
end
    


