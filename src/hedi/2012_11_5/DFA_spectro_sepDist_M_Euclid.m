%% Cette routine effectue une PCA et une DFA sur les signaux de zebra finch (distance calls) à 2m. 
% Changer le path pour faire la même analyse aux autres distances : analyse distances séparées
% Seuls 15 cris sont utilisés, les 16èmes cris pour chaque oiseau sont le validating dataset pour la cross-validation.

%% Spectrogrammes de cris à 2m pour récupérer les param. les plus importants
% Obtention des données de spectros à 2m

clear all
close all
ndist = input('Distance de propagation (en m) ?: ');
cd('/Users/Solveig/PhD/PROPAG/PROPAG BASE DE DONNEES/Selection Propagmale2m')

fs=44100;       % Sampling frequency
fband = 70;    % Frequency band for the spectrogram A DIMINUER pour augmenter résolution
dBscale = 50;   % Decibel scale for thresholding the spectrogram
nbirds=16;      % Total number of birds
ncalls=16;      % Total number of calls per bird

fmax = 8000;   % fréquence maximale pour le filtre (corresp. à l'audiogramme du Z.f)
fmin = 500;    % fréquence minimale pour le filtre

nb = 25;        % nombre de PC considérés A VOIR SI ON MODIFIE !
nDF = 4;        % nombre de DF considérées

coordcri = load('coordcriM0m.csv'); 
MaxLength = 0.29090909; % length of the longest call, mbird2call3

isound1 = wavread('mbird2call3.wav'); % on charge le cri le plus long
debcri = coordcri(19,3); % coordonnées de Mbird2call3, début
fincri = coordcri(19,4); % coordonnées de Mbird2call3, fin
% isound1_cut = isound1(fix(debcri*fs-0.001*fs):fix(fincri*fs+0.001*fs));
isound1_cut = isound1(fix(debcri*fs):fix(fincri*fs));
maxlen = length(isound1_cut); % en sampling points

figure(1);
[t1, f1, logS1] = specSlvg(isound1_cut, fband, fs, dBscale);
indMax = find(f1>=fmax, 1); % première ligne de f1 avec valeur > fmax
indMin = find(f1<=fmin, 1);
logS1 = logS1(indMin:indMax, :); % on réduit logS1 aux bornes de filtration

[nf, nt] = size(logS1);

X=zeros(nbirds*ncalls, nt*nf); % créer matrice totale de résultats
Group=zeros(nbirds*ncalls, 1); % créer matrice des noms d'individus

dirName = sprintf('/Users/Solveig/PhD/PROPAG/PROPAG BASE DE DONNEES/Selection Propagmale%dm', ndist);
cd(dirName)
for ib = 1:nbirds
    for ic = 1:ncalls
        fileName = sprintf('mbird%dcall%d.wav', ib, ic);
        icall = (ib-1)*nbirds + ic;
        isound_full = wavread(fileName);
        %fprintf('%s\n', fileName);
        
        % On ne prend que le cri et on le centre dans une fenêtre égale au
        % cri le plus long -> même longeur de fichier pour tt le monde
        clear debcri fincri;
        debcri = coordcri(16*(ib-1)+ic,3);
        fincri = coordcri(16*(ib-1)+ic,4);
        isound_cut = isound_full(fix(debcri*fs):fix(fincri*fs));
        isound = zeros(maxlen, 1); % taille du fichier égale à celle du son le plus long
        decal = length(isound)-length(isound_cut); % différence avec son le plus long
        if decal <0
            fprintf('Cri plus long que durée maximale, impossible !\n');
            pause();
        end
        if decal == 0
            isound = isound_cut; % cas particulier du cri le plus long
        else
            % isound(fix(decal/2):fix(decal/2)+length(isound_cut)-1)=isound_cut;
            % au dessus : avec silence de part et d'autre du cri
            isound=isound_full(fix(debcri*fs-decal/2):fix(fincri*fs+decal/2));
        % on place le signal coupé au "milieu" du vecteur à taille du cri
        % le plus long, en gardant le bruit de fond
        end
        
        % plot the spectrogram
        figure(2);
        clear t f logS;
        [t, f, logS] = specSlvg(isound, fband, fs, dBscale); % !! fonction spec() modifiée pour 
        % replacer les -Inf par le min value au dessus de dBscale
        logS = logS(indMin:indMax, :);
        
        %badind = find(logS == -Inf);
        
        if size(logS) ~= [nf nt]
            fprintf('Unexpected size for file "mbird%dcall%d"\n', ib, ic);
            return;
        end % vérifier que tous les fichiers ont la même taille
        
        X(icall, :) = reshape(logS, 1, nt*nf);
        Group(icall)=ib;
    %pause(); % pause pour tous les cris
    end
    %pause(); % pause entre chaque bird
end

Xmoy = mean(X);

%save('X_data', 'X');
%save('Group_data', 'Group');
%save('X_mean', 'Xmoy');

% On modifie le vecteur Group pour enlever une ligne par oiseau (ATTENTION ne marche que pour fitting dataset = 15)
% Le nouveau vecteur Group n'aura que 15 lignes par oiseau, pour les 15 cris sélectionnés comme fitting dataset.

Group_fit = Group;
for i = 1:nbirds % on enlève la première ligne de chaque oiseau
    if i ==1
        Group_fit(i) = [];
    else
        Group_fit((i-1)*nbirds) = [];
    end
end

%% Début boucle ! Sélection du fitting dataset et validating dataset
% On ne prend que 15 cris par oiseau, au hasard, pour le fitting dataset

it_num = 250; % nombre d'itérations à changer éventuellement
ConfMat = zeros(nbirds, nbirds); % matrice de confusion totale
PCC_M = zeros(it_num, 1); % PCC moyen calculé à chaque itération
PCCbird = zeros(nbirds,1); % Classif correcte pour chaque oiseau, permet de calculer PCC

for it = 1:it_num 
    
    % Matrice avec numéro du cri pour chaque oiseau qui sera utilisé pour validating dataset
    X_fit = X;
    CallOut = randi(ncalls, 1, nbirds);
    CallInd = zeros(1, nbirds);

    for i = 1:nbirds % on remplace les lignes corresp aux validating calls par NaN
        X_fit((i-1)*nbirds + CallOut(i),:) = NaN;
        CallInd(i) = (i-1)*nbirds + CallOut(i);
    end

    for i = 1:nbirds % on enlève les lignes corresp. aux validating calls
        if i ==1
            X_fit(CallInd(i),:) = [];
        else
            X_fit(CallInd(i)-i+1,:) = [];
        end
    end

    for i = 1:nbirds*ncalls-nbirds % On vérifie l'intégrité de X (pas de NaN)
            NaN_ind = find(isnan(X_fit));
            if length(NaN_ind) >0
            fprintf('\nSTOP: %d NaN values found in fitting dataset (X_fit), there is a glitch!\n', length(NaN_ind));
            return
            end
    end

    [nr1 nc1] = size(X_fit); % vérification de compatibilité entre Group_fit et X_fit
    [nr2 nc2] = size(Group_fit);
    if nr1 ~= nr2
        fprintf('\n STOP: Group and fitting dataset (X_fit) do not have the same size!\n');
        return
    end

    [Coeff, Score, latent] = princomp(X_fit); % PCA pour le fitting dataset
    [d, p, stats] = manova1(Score(:, 1:nb), Group_fit); % Manova sur la fitting dataset

    [mean_grp, sem_grp, meanCI_grp, range_grp] = grpstats(stats.canon(:,1:nDF),Group_fit, {'mean', 'sem', 'meanci', 'range'});

    %fprintf('Eigen values pour chaque DF :\n');
    %disp(stats.eigenval);
    %fprintf('On considère les %d 1eres DF pour la suite. Si not OK, changer nDF.\n', nDF)

    % Application du validating dataset aux DF obtenues à 2m
    % Création d'une matrice avec les cris sélectionnés comme validating dataset
    for i = 1:nbirds
        X_valid(i,:) = X((i-1)*nbirds + CallOut(i),:);
    end

    if size(X_valid) ~= [nbirds nt*nf] % Vérification taille X_valid
        fprintf('\n STOP: validating dataset has a wrong size!\n');
        return
    end

    % Application du validating dataset au fitting dataset
    X_valid_moy = mean(X_valid);
    Xc = X_valid - repmat(X_valid_moy, size(X_valid,1), 1); % matrice X centrée
    Y = Xc*Coeff(:, 1:nb); % Y correspond au Score pour le validating dataset (size = 256, 20)
    Ymoy = mean(Y);
    Yc = Y - repmat(Ymoy, size(Y,1), 1); % matrice Y centrée
    Canon = Yc*stats.eigenvec(:, 1:nDF); % Canon à la distance voulue pour les 3 premières DF (si nDF=3)

    % Calcul des distances Euclidiennes
    Mah = zeros(nbirds, nbirds); % ATTENTION ne marche que pour validating dataset de 16 cris !!
    for i = 1:nbirds
        for j = 1:nbirds
            Mah(i,j) = sqrt((Canon(i,:) - mean_grp(j,:))*(Canon(i,:) - mean_grp(j,:))');
        end
    end

    % Confusion Matrix : addition des classif à chaque itération
    for i = 1:nbirds
            clear row col minVal
            if min(Mah(i,:)) == Mah(i, i)
                ConfMat(i,i) = ConfMat(i,i)+1; % classification correcte
            else
                minVal = min(Mah(i,:));
                [row, col] = find(Mah(i,:)==minVal);
                ConfMat(i,col) = ConfMat(i,col)+1; % classification incorrecte
            end
    end
    for i = 1:nbirds
    PCCbird(i) = ConfMat(i,i) / sum(ConfMat(i, :), 2) * 100; % sum(.., 2) = somme par ligne
    end
    PCC_M(it) = mean(PCCbird);
    fprintf('Iteration #%d\n', it);
end

%% Wrap up figure et sauvegarde résultats

% Confusion Matrix
figure();
imagesc(ConfMat)
colormap(gray)
colorbar
title(sprintf('Confusion matrix %dm %.1f%% Correct', ndist, 100*sum(diag(ConfMat))/(nbirds*it_num)));

% Calcul des pourcentages de classification correcte
CorrClass = zeros(nbirds, 2);

for i = 1:nbirds
    CorrClass(i, 1) = i;
    CorrClass(i, 2) = ConfMat(i,i) / sum(ConfMat(i, :), 2) * 100; % sum(.., 2) = somme par ligne
end

PCCmoy = mean(CorrClass(:,2)); % PCC moyens sur tous les oiseaux
PCCsd = std(CorrClass(:,2));

fprintf('\nPourcentage moyen de classification correcte à %dm : %.2f%%\n', ndist, PCCmoy);
fprintf('Ecart-type des PCC : %.2f\n\n', PCCsd);

DataSpectroDFA = struct;
DataSpectroDFA.Dist = sprintf('%dm', ndist);
DataSpectroDFA.ConfMat = ConfMat;
DataSpectroDFA.CorrClass = CorrClass;
DataSpectroDFA.PCCit = PCC_M;

% Sauvegarde de la structure de données sur le disque dur
cd('/Users/Solveig/PhD/PROPAG/aRESULTS PROPAG/IMAGES DFA SPECTROS/SPECTROS_sep_dist_M/');
save('DataSpectroDFA2mMeuclid', 'DataSpectroDFA')

%% Représentation des DF1 à 4 sous forme de spectros : ce que les axes veulent dire
% Ce calcul est indépendant de la boucle de cross-validation.

% On refait l'analyse sur toute la matrice X (tous cris confondus) :
[CoeffTot, ScoreTot, latentTot] = princomp(X); % PCA pour tout le dataset
[dTot, pTot, statsTot] = manova1(ScoreTot(:, 1:nb), Group); % Manova (Canon is inside 'stats' structure)

PC_DF = CoeffTot(:, 1:nb) * statsTot.eigenvec(:, 1:nDF); % les quatre 1ères DF se lisent en colonnes

% définir les min et max des futures DF pour donner les mêmes limites de
% colorscale à tous les graphiques issus de imagesc()
clear cmin cmax clims
cmin = min(min(PC_DF(:, 1:nDF)));
cmax = max(max(PC_DF(:, 1:nDF)));
cabs = max(abs(cmin), abs(cmax));
clims = [-cabs cabs];

figure();
PC_DF_1 = reshape(PC_DF(:,1), nf, nt);
imagesc(t, f(indMin:indMax), PC_DF_1, clims)
title(sprintf('DF1 (using %d PC)', nb));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');

figure();
PC_DF_2 = reshape(PC_DF(:,2), nf, nt);
imagesc(t, f(indMin:indMax), PC_DF_2, clims)
title(sprintf('DF2 (using %d PC)', nb));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');

figure();
PC_DF_3 = reshape(PC_DF(:,3), nf, nt);
imagesc(t, f(indMin:indMax), PC_DF_3, clims)
title(sprintf('DF3 (using %d PC)', nb));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');

figure();
PC_DF_4 = reshape(PC_DF(:,4), nf, nt);
imagesc(t, f(indMin:indMax), PC_DF_4, clims)
title(sprintf('DF4 (using %d PC)', nb));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
    

