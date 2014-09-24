function [calfilename] = RandomPredictConfusionMatrices_SM(MatfilePath, Birdname)

Bootmax = 1000;      %Number of random matrices to compute per cell
ndist = 5;
load(MatfilePath);  
%% Retrieve the best Individual confusion matrix of that cell
% tu as besoin de faire cette étape (retrouver la matrice avec la fenêtre que tu as choisi) uniquement pour les matrices calculées
% à 2m puisque pour toutes les autres distances tu ne t'es concentrée que sur une seule fenêtre (fenêtre de calcul de la matrice de confusion).
% Ensuite tous le calcul dans ce code concerne le bootstrap pour une matrice de confusion par cellule, toi tu as plusieurs distances + 2
% catégories de cris (syn vs calls) donc il faut que tu fasse un loop pour faire tourner ce code sur toutes les matrices.
% SM : pour l'instant je me focalise sur les data Con et à 2m (avec la bonne winSize)

% Take care of 2m first (got get the max window size)
wind = data{1}{1}(1).ChosenWinSize;
index = find(data{1}{1}(1).winsize==wind);
VocConfMat = data{1}{1}(index).ConfMatTot;
VocTypeSel = data{1}{1}(index).StimID_2m(:,1);
birdID=unique(VocTypeSel);
NbirdID=length(birdID);
% ncalls = length(VocTypeSel)/NbirdID;
% if (ncalls ~= round(ncalls))
%     fprintf('Error in calculating the number of calls for each individual: the number is not an integer!\n');
%     return;
% end

ConfMat_Rand = cell(Bootmax, 1);
List_VocRand = cell(Bootmax, 1);
ConfMat_Bird = cell(Bootmax, 1);
infoRand_perBird = zeros(Bootmax, 1);

for bb=1:Bootmax
    %fprintf('%d\n', bb);
    rng('shuffle');
    VocTypeSel_rand=VocTypeSel(randperm(length(VocTypeSel)));
    RR=0;
    VocRand = zeros(1,length(VocTypeSel));
    NVPC = zeros(NbirdID,1);
    for vtR=1:NbirdID
        stR=birdID(vtR);
        selectorR=strcmp(VocTypeSel_rand, stR);
        selectorIndR=find(selectorR);
        VocRand(1,RR+1:RR+length(selectorIndR)) = selectorIndR;
        RR=RR+length(selectorIndR);
        NVPC(vtR) = length(selectorIndR);
    end
    confusion_matrix_vocalizationsRand = VocConfMat(:, VocRand);
    confusion_matrix_vocalizationsRand = confusion_matrix_vocalizationsRand./sum(sum(confusion_matrix_vocalizationsRand)); % cette ligne n'est pas nécessaire si tu travaille avec la matrice de proba jointes
    ConfMat_Rand{bb} = confusion_matrix_vocalizationsRand; % cette cellule sauve toutes les matrices random
    List_VocRand{bb} = VocRand;
    
    ncalls = NVPC(1);
    ConfMatTemp=zeros(NbirdID);
    clear BirdMat
    for row=1:NbirdID
        for col=1:NbirdID
            BirdMat=confusion_matrix_vocalizationsRand((row-1)*ncalls+1:(row-1)*ncalls+ncalls, (col-1)*ncalls+1:(col-1)*ncalls+ncalls);
            ConfMatTemp(row,col) = sum(sum(BirdMat));
        end
    end
    if (sum(sum(ConfMatTemp))<0.9999 || sum(sum(ConfMatTemp))>1.0001) % the sum is not always completely equal to 1 for some reason (= to "1.0000"...)
        error('ERROR: The sum of probabilities in the confusion matrix is not equal to 1 in ConfMatTemp!');
    end
    mi_confTemp = info_matrix(ConfMatTemp);
    infoRand_perBird(bb) = mi_confTemp;  % mutual information for the Confusion Matrix calculated per bird (per individual to discriminate)
    ConfMat_Bird{bb} = ConfMatTemp;      % Confusion Matrix calculated per bird
end

RandMat.subject = Birdname;
RandMat.originalfile = MatfilePath;
RandMat.mi_perBird{1} = infoRand_perBird; % for 2 m
RandMat.ConfMat_Bird{1} = ConfMat_Bird;
RandMat.ConfMat_Rand{1} = ConfMat_Rand;
RandMat.List_VocRand{1} = List_VocRand;

%% Loop through all other distances

for di = 2:ndist
    VocConfMat = data{1}{di}.ConfMatTot;
    
    ConfMat_Rand = cell(Bootmax, 1);
    List_VocRand = cell(Bootmax, 1);
    ConfMat_Bird = cell(Bootmax, 1);
    infoRand_perBird = zeros(Bootmax, 1);
    
    for bb=1:Bootmax
        rng('shuffle');
        VocTypeSel_rand=VocTypeSel(randperm(length(VocTypeSel)));
        RR=0;
        VocRand = zeros(1,length(VocTypeSel));
        NVPC = zeros(NbirdID,1);
        for vtR=1:NbirdID
            stR=birdID(vtR);
            selectorR=strcmp(VocTypeSel_rand, stR);
            selectorIndR=find(selectorR);
            VocRand(1,RR+1:RR+length(selectorIndR)) = selectorIndR;
            RR=RR+length(selectorIndR);
            NVPC(vtR) = length(selectorIndR);
        end
        confusion_matrix_vocalizationsRand = VocConfMat(:, VocRand);
        confusion_matrix_vocalizationsRand = confusion_matrix_vocalizationsRand./sum(sum(confusion_matrix_vocalizationsRand)); % cette ligne n'est pas nécessaire si tu travaille avec la matrice de proba jointes
        ConfMat_Rand{bb} = confusion_matrix_vocalizationsRand; % cette cellule sauve toutes les matrices random
        List_VocRand{bb} = VocRand;
        
        ncalls = NVPC(1);
        ConfMatTemp=zeros(NbirdID);
        clear BirdMat
        for row=1:NbirdID
            for col=1:NbirdID
                BirdMat=confusion_matrix_vocalizationsRand((row-1)*ncalls+1:(row-1)*ncalls+ncalls, (col-1)*ncalls+1:(col-1)*ncalls+ncalls);
                ConfMatTemp(row,col) = sum(sum(BirdMat));
            end
        end
        if (sum(sum(ConfMatTemp))<0.9999 || sum(sum(ConfMatTemp))>1.0001) % the sum is not always completely equal to 1 for some reason (= to "1.0000"...)
            error('ERROR: The sum of probabilities in the confusion matrix is not equal to 1 in ConfMatTemp!');
        end
        mi_confTemp = info_matrix(ConfMatTemp);
        infoRand_perBird(bb) = mi_confTemp;  % mutual information for the Confusion Matrix calculated per bird (per individual to discriminate)
        ConfMat_Bird{bb} = ConfMatTemp;      % Confusion Matrix calculated per bird
    end
    
    RandMat.mi_perBird{di} = infoRand_perBird; % for 2 =
    RandMat.ConfMat_Bird{di} = ConfMat_Bird;
    RandMat.ConfMat_Rand{di} = ConfMat_Rand;
    RandMat.List_VocRand{di} = List_VocRand;
end

[Path, Matfile] = fileparts(MatfilePath);
if ismac()
        [status username] = system('who am i');
        if strcmp(strtok(username), 'Solveig')
            calfilename = fullfile('/Users','Solveig','PhD','ELECTROPHY','Neuro_DATA','RandMat_bootstrap',['RandMat_' Matfile '.mat']);
        end
else
    calfilename = fullfile('/auto','k8','fdata', 'solveig','matfile',['RandMat_' Matfile '.mat']);
end
save(calfilename, 'RandMat');
fprintf(1,'\nDone making calculus on %s\nData saved under %s\n', MatfilePath, calfilename);
end



