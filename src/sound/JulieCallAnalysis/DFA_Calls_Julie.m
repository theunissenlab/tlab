%% Perform a regularized DFA based on Spectrogrammes for all the zebra finch vocalizations

% Read the data base produced by VocSectioningAmp.m 
load /Users/frederictheunissen/Documents/Data/Julie/BanqueCrisNeuro/vocCuts.mat

% To make RAM space, clear the sound files
clear soundCutsTot;

% chose to normalize or not...
normFlg = 1;
nSample = size(spectroCutsTot, 1);

if normFlg % Normalize by peak
    for i=1:nSample
        spectroCutsTot(i,:) = spectroCutsTot(i,:) - max(spectroCutsTot(i,:));
    end
end

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
if normFlg
    save /Users/frederictheunissen/Documents/Data/Julie/BanqueCrisNeuro/vocPCANorm.mat
else
    save /Users/frederictheunissen/Documents/Data/Julie/BanqueCrisNeuro/vocPCA.mat
end

%% Start here if you already done the PCA
%load /Users/frederictheunissen/Documents/Data/Julie/BanqueCrisNeuro/vocPCA.mat
load /Users/frederictheunissen/Documents/Data/Julie/BanqueCrisNeuro/vocPCANorm.mat

%% Display the first nb PCs as Spectrograms
figure(1);
nb = 100; % Show the first 100 PCs here.

clear cmin cmax clims
cmin = min(min(Coeff(:, 1:nb)));
cmax = max(max(Coeff(:, 1:nb)));
cabs = max(abs(cmin), abs(cmax));
clims = [-cabs cabs];
nf = length(fo);
nt = length(to);


for i=1:nb
    subplot(ceil(nb/10),10,i);
    PC_Spect = reshape(Coeff(:,i), nf, nt);
    imagesc(to, fo, PC_Spect, clims)
    title(sprintf('PCA%d', i));
    axis xy;
    axis off;
end

figure(2);
nb = 300;  % Show the cumulative distribution up to 300 pcs.
cumlatent = cumsum(latent);
plot(cumlatent(1:nb)./sum(latent));
xlabel('PCA Dimensions');
ylabel('Cumulative Variance');


%%  Perform the classical DFA on all of the Data to obtain DFA for vocalization Type


% A sanity check using randomization:
% nr = length(vocTypeCuts);
% [nDF, p, stats] = manova1(Score(:, 1:nb), vocTypeCuts(randperm(nr)));
% This indeed gave 0 significant dimensions

% Calculate first for 300 PCs
nb = 300;
[nDF300, p, stats] = manova1(Score(:, 1:nb), vocTypeCuts);
[mean_grp, sem_grp, meanCI_grp, range_grp, name_grp] = grpstats(stats.canon(:,1:nDF300),vocTypeCuts', {'mean', 'sem', 'meanci', 'range', 'gname'});
ngroups = size(mean_grp,1);
PC_DF_300 = Coeff(:, 1:nb) * stats.eigenvec(:, 1:nDF300);

DF_corr = zeros(5, nDF300);

for nbi=1:5
    nb = nbi*50;   % number of PCs
    [nDF, p, stats] = manova1(Score(:, 1:nb), vocTypeCuts);
    [mean_grp, sem_grp, meanCI_grp, range_grp, name_grp] = grpstats(stats.canon(:,1:nDF),vocTypeCuts', {'mean', 'sem', 'meanci', 'range', 'gname'});
    ngroups = size(mean_grp,1);
    
    %  Display the significant DFA in the spectrogram space
    PC_DF = Coeff(:, 1:nb) * stats.eigenvec(:, 1:nDF);
    
    for i = 1:nDF300
            if i <= nDF
                DF_corr(nbi,i) = corr(PC_DF(:,i), PC_DF_300(:,i));
            end
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
    nDFPlot = 10;
    for i=1:nDF
        if (i <= nDFPlot)
            subplot(6,nDFPlot,i+(nbi-1)*nDFPlot);
            PC_Spect = reshape(PC_DF(:,i), nf, nt);
            imagesc(to, fo, PC_Spect, clims)
            if i == 1
                title(sprintf('Voc Type DFs using %d PC)', nb));
            end
            axis xy;
            axis([0 0.2 0 8000]);
            axis off;
        end
    end
end

% Plot the last line
for i=1:nDF300
    if (i <= nDFPlot)
        subplot(6,nDFPlot,i+5*nDFPlot);
        PC_Spect = reshape(PC_DF_300(:,i), nf, nt);
        imagesc(to, fo, PC_Spect, clims)
        if i == 1
            title(sprintf('Voc Type DFs using %d PC)', 300));
        end
        axis xy;
        axis([0 0.2 0 8000]);
        axis off;
    end
end

% Plot the coefficients of the DFA for n=50, n=100, n=200 PCs
figure(4);
% Using 50 PCs
subplot(1,4,1);
[nDF, p, stats] = manova1(Score(:, 1:50), vocTypeCuts); 
weigthsPC_DFA = zeros(300, nDF);
weigthsPC_DFA(1:50, :) = abs(stats.eigenvec(1:50, 1:nDF));
imagesc(1:nDF, 1:300, weigthsPC_DFA);
ylabel('PCA Vector');
title('50 PCs');

% Using 100 PCs
subplot(1,4,2);
[nDF, p, stats] = manova1(Score(:, 1:100), vocTypeCuts); 
weigthsPC_DFA = zeros(300, nDF);
weigthsPC_DFA(1:100, :) = abs(stats.eigenvec(1:100, 1:nDF));
imagesc(1:nDF, 1:300, weigthsPC_DFA);
ylabel('PCA Vector');
title('100 PCs');

% Using 200 PCs
subplot(1,4,3);
[nDF, p, stats] = manova1(Score(:, 1:200), vocTypeCuts); 
weigthsPC_DFA = zeros(300, nDF);
weigthsPC_DFA(1:200, :) = abs(stats.eigenvec(1:200, 1:nDF));
imagesc(1:nDF, 1:300, weigthsPC_DFA);
ylabel('PCA Vector');
title('200 PCs');

% Using 300 PCs
subplot(1,4,4);
[nDF, p, stats] = manova1(Score(:, 1:300), vocTypeCuts); 
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
% THe problem with this bootstrap is that measures are not indenpendent
% since many samples of each come from same bird

n_boot = 1000;     % Number of permutations
n_valid = fix(0.1*ncutsTot); % Save 10% of the data for cross-validation

PCC_info = struct('nb', {50, 100, 150, 200, 250, 300}, 'PCC_Total', 0, 'PCC_M', 0, 'PCC_group', zeros(ngroups, ngroups));

for inb = 1:length(PCC_info)
    
    % Number of PCA used in the DFA
    nb = PCC_info(inb).nb;
    
    % Allocate space for distance vector and confusion matrix
    Dist = zeros(1, ngroups);
    ConfMat = zeros(ngroups, ngroups);
    
    for iboot=1:n_boot
        
        ind_valid = randsample(ncutsTot, n_valid);    % index of the validation calls
        
        % Separate data into fitting and validation
        X_valid = Score(ind_valid, 1:nb);
        X_fit = Score(:, 1:nb);
        X_fit(ind_valid, :) = [];
        
        % Similarly for the group labels.
        Group_valid = vocTypeCuts(ind_valid);
        Group_fit = vocTypeCuts;
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
legend('Total', 'Averaged per Type');



%%  Repeat the boot-strap by systematically swapping birds

birdNames = unique(birdNameCuts);
nBirds = length(birdNames);
vocTypes = unique(vocTypeCuts);   % This returns alphabetical but it is not the same as the order returned by grpstats above
ngroups = length(vocTypes);
nPerm = 1000;



PCC_info_bird = struct('nb', {10, 20, 30, 40, 50, 75, 100, 150, 200, 250, 300}, 'PCC_Total', 0, 'PCC_M', 0, 'PCC_group', zeros(ngroups, ngroups));

for inb = 1:length(PCC_info_bird)
    
    % Number of PCA used in the DFA
    nb = PCC_info_bird(inb).nb;
    
    % Allocate space for distance vector and confusion matrix
    Dist = zeros(1, ngroups);
    ConfMat = zeros(ngroups, ngroups);
    n_validTot = 0;
    
    for iperm=1:nPerm
        
        birdForEachVoc = randperm(nBirds, ngroups);
        ind_valid = [];
        for igroup = 1:ngroups
            ind_valid = [ind_valid find(strcmp(birdNameCuts, birdNames{birdForEachVoc(igroup)}) & strcmp(vocTypeCuts, vocTypes{igroup})) ];
        end
        
        % ind_valid = find(strcmp(birdNameCuts, birdNames{ibird}));    % index of the validation calls
        n_valid = length(ind_valid);
        
        % Separate data into fitting and validation
        X_valid = Score(ind_valid, 1:nb);
        X_fit = Score(:, 1:nb);
        X_fit(ind_valid, :) = [];
        
        % Similarly for the group labels.
        Group_valid = vocTypeCuts(ind_valid);
        Group_fit = vocTypeCuts;
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
        
        n_validTot = n_validTot + n_valid;
        
    end
    
    PCC_Total = 100.0*sum(diag(ConfMat))./n_validTot;
    PCC_group = zeros(ngroups, ngroups);
    for i = 1:ngroups
        for j = 1:ngroups
            PCC_group(i,j) = ConfMat(i,j) / sum(ConfMat(i, :), 2) * 100; % sum(.., 2) = somme par ligne
        end
    end
    PCC_M = mean(diag(PCC_group));
    
    % Store the information
    
    PCC_info_bird(inb).PCC_Total = PCC_Total;
    PCC_info_bird(inb).PCC_M = PCC_M;
    PCC_info_bird(inb).PCC_group = PCC_group;
    
end
% Display confusion Matrix

figure(6);
for inb = 1:length(PCC_info_bird)
    subplot(1,length(PCC_info_bird),inb);
    imagesc(PCC_info_bird(inb).PCC_group);
    xlabel('Guess');
    ylabel('Actual');
    colormap(gray);
    colorbar;
    title(sprintf('Confusion Matrix RF %.1f%%(%.1f%%) Correct', PCC_info_bird(inb).PCC_Total, PCC_info_bird(inb).PCC_M));
    set(gca(), 'Ytick', 1:ngroups);
    set(gca(), 'YTickLabel', name_grp);
    set(gca(), 'Xtick', 1:ngroups);
    set(gca(), 'XTickLabel', name_grp);
end

figure(7);
plot([PCC_info_bird.nb], [PCC_info_bird.PCC_Total],'k');
hold on;
plot([PCC_info_bird.nb], [PCC_info_bird.PCC_M],'r');
xlabel('Number of PCs');
ylabel('Correct Classification');
legend('Total', 'Averaged per Type');
hold off;

% save the DFA
if normFlg
    save vocTypeDFANormBird.mat fo to PC_DF PCC_info_bird ngroups name_grp
else
    save vocTypeDFABird.mat fo to PC_DF PCC_info_bird ngroups name_grp
end

%% From these results we choose nb=50 PCs and make nice plots + save some data
% You can load vocTypeDFA.mat and vocPCA.mat or the normalized versions and and start here to make these
% final plots
%load /Users/frederictheunissen/Documents/Data/Julie/BanqueCrisNeuro/vocPCANorm.mat

%load /Users/frederictheunissen/Documents/Data/Julie/analysis/vocTypeDFANorm.mat

inb = 5;       % Corresponds to 50 PCs in PCC_info_bird
nb = 50;


% A nice confusion matrix.
figure(8);

imagesc(PCC_info_bird(inb).PCC_group);
xlabel('Guess');
ylabel('Actual');
colormap(gray);
colorbar;
title(sprintf('nPCA = %d PC Total= %.1f%% PCC Mean = %.1f%%', PCC_info_bird(inb).nb, PCC_info_bird(inb).PCC_Total, PCC_info_bird(inb).PCC_M));
set(gca(), 'Ytick', 1:ngroups);
set(gca(), 'YTickLabel', name_grp);
set(gca(), 'Xtick', 1:ngroups);
set(gca(), 'XTickLabel', name_grp);

figure(9);

[nDF, p, statsDFA] = manova1(Score(:, 1:nb), vocTypeCuts);
[mean_grp, std_grp, name_grp] = grpstats(statsDFA.canon(:,1:nDF),vocTypeCuts', {'mean', 'std', 'gname'});

    
%  Display the significant DFA in the spectrogram space
PC_DF = Coeff(:, 1:nb) * statsDFA.eigenvec(:, 1:nDF);

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
    
    title(sprintf('DFA %d', i));
    
    axis xy;
    axis([0 0.2 0 8000]);
    if i ~= 1
        axis off;
    end
    
    if (i == 1)
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
    end
end

% Make DF1 vs DFX plots
figure(10);
name_grp_plot = {'Be', 'LT', 'Tu', 'Th', 'Di', 'Ag', 'Wh', 'Ne', 'Te', 'DC', 'So'};
    
colorVals = { [0 230 255], [0 95 255], [255 200 65], [255 160 60], [255 105 15],...
    [255 0 0], [255 0 255], [255 100 255], [255 180 255], [140 100 185], [0 0 0]}; 

for iDF=2:nDF
    subplot(2,ceil((nDF-1)/2), iDF-1);
    for ig=1:ngroups
        for ig_ind=1:ngroups
            if strcmp(name_grp_plot{ig}, name_grp{ig_ind})
                break;
            end
        end
        
        plot(mean_grp(ig_ind,1), mean_grp(ig_ind,iDF), 's', 'MarkerSize', 10, ...
                'color', colorVals{ig}./255,'MarkerEdgeColor','k',...
                'MarkerFaceColor',colorVals{ig}./255);
        hold on;
    end
    xlabel('DF1');
    ylabel(sprintf('DF%d', iDF));
    axis([-5 5 -5 5]);
    axis square;
    if iDF == nDF
        legend(name_grp_plot);
    end
    hold off;
end
    


%% Repeat the classification with the Random Forest Algorithm   

% The random tree does its own bootstrat and we don't need a training and
% validation set.  
nb = 50;                    % Here we choose 50 PCs as well to match DFA

Group_names = unique(vocTypeCuts);
ngroups = length(Group_names);

ConfMat = zeros(ngroups, ngroups);   

% THis worked on a prior version of Matlab... B = TreeBagger(300, Score(:, 1:nb), vocTypeCuts, 'OOBPred', 'on', 'priorprob', 'equal', 'MinLeaf', 5, 'NPrint', 10);
B = TreeBagger(300, Score(:, 1:nb), vocTypeCuts, 'FBoot', 1.0, 'OOBPred', 'on', 'MinLeaf', 5, 'NPrint', 10);

Group_predict = oobPredict(B);   % This returns the predictions for the out of bag values.

n_valid = length(vocTypeCuts);   % this is the total number of observations

for i = 1:n_valid
    k_actual = find(strcmp(Group_names,vocTypeCuts(i)));
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
set(gca(), 'YTickLabel', Group_names);
set(gca(), 'Xtick', 1:ngroups);
set(gca(), 'XTickLabel', Group_names);

%% Perform the same analysis using multinomial regression


% Nultinomial regression is slow and only converges for a small number of
% PCs (should be less than minimum samples in each group)
% Maybe we should do this with the DFA coefficient? (stats.canon?)

Group_names = unique(vocTypeCuts);
ngroups = length(Group_names);

nScores = size(Score,1);
yVocType = zeros(nScores, 1);
nSampleGroup = zeros(ngroups, 1);

for ig = 1:ngroups
    ind = find(strcmp(vocTypeCuts, Group_names{ig}));
    yVocType(ind) = ig;
    nSampleGroup(ig) = length(ind);
    
end

% Fitting the multinomial model
nb = (min(nSampleGroup) - 1)./2;
[BM, devianceM, statsM] =  mnrfit(Score(:, 1:nb), yVocType);

PC_MN = Coeff(:, 1:nb) * BM(2:end, :);

%% Perform the logistic regression in a pairwise fashion
% Here we used the space obtained in the DFA...

Group_names = unique(vocTypeCuts);
ngroups = length(Group_names);
nScores = size(Score,1);

nb = 50;     % should match the nb used to obtain the statsDFA used below
PC_LR = zeros(size(Coeff,1),ngroups);

for ig = 1:ngroups
    Group_names{ig}
    yVocType = zeros(nScores, 1);
    yVocType(strcmp(vocTypeCuts, Group_names{ig})) = 1;
    
    [BLR, devianceLR, statsLR] =  glmfit(statsDFA.canon(:,1:nDF), yVocType, 'binomial');
    PC_LR(:,ig) = Coeff(:, 1:nb) * (statsDFA.eigenvec(:,1:nDF)*BLR(2:end, :));
    
end

figure(13);
nf = length(fo);
nt = length(to);
for ig = 1:ngroups
    clear cmin cmax clims
    cmin = min(PC_LR(:,ig));
    cmax = max(PC_LR(:,ig));
    cabs = max(abs(cmin), abs(cmax));
    clims = [-cabs cabs];
    subplot(1, ngroups, ig);
    imagesc(to,fo, reshape(PC_LR(:, ig), nf, nt), clims);
    axis xy;
    axis off;
    title(Group_names{ig});
end
