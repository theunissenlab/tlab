cd /auto/k6/julie/matfile
% resultsDirectory='/auto/k6/julie/matfile';
% addpath('/auto/k1/queued');


input_dir=pwd;
Subjects = dir(input_dir);
%cmd = 'calls_selectivity_ConfusionMatrix(''%s'');';
optWin=zeros(1410,1);
optWin_BGwo=zeros(1410,1);
MI_confB=zeros(1410,9);
 MI_confBBGwo=zeros(1410,9);
ii=0;

for ss=1:length(Subjects)
    Indiv=Subjects(ss).name;
    if length(Indiv)==11
%     if strcmp(Indiv, 'GreBlu9508M')
        Idir=fullfile(input_dir, Indiv);
        matfiles=dir(fullfile(Idir,'ConfMat*.mat'));
%        matfiles=dir(fullfile(Idir,'FirstVoc*Site2*ss1*.mat'));
        lm=length(matfiles);
        SS_Ind=zeros(lm,1);
        for ff=1:lm
            if ~isempty(strfind(matfiles(ff).name, 'ss'))
                SS_Ind(ff)=1;
            end
            %if ~isempty(strfind(matfiles(ff).name, 'e12'))
%                 SS_Ind(ff)=1;
%             end
%             if ~isempty(strfind(matfiles(ff).name, 'e21'))
%                 SS_Ind(ff)=1;
%             end
        end
        Indices=find(SS_Ind);
        LM=length(Indices);
        for hh=1:LM
            
            MatfilePath=fullfile(Idir, matfiles(Indices(hh)).name)
            if isempty(strfind(matfiles(Indices(hh)).name, 'ConfMat_Site2_L1100R1450_e12_s0_ss')) && isempty(strfind(matfiles(Indices(hh)).name, 'ConfMat_Site2_L1100R1450_e21_s0_ss'))
                ii=ii+1;
                MAT=load(MatfilePath);
                Winsize=MAT.winSize;
                MAXWin=find(MAT.mi_confusionB==max(MAT.mi_confusionB));
                MAXWin_BGwo=find(MAT.mi_confusionB_BGwo==max(MAT.mi_confusionB_BGwo));
                optWin(ii)=Winsize(MAXWin);
                optWin_BGwo(ii)=Winsize(MAXWin_BGwo);
                MI_confB(ii,:)=MAT.mi_confusionB;
                MI_confBBGwo(ii,:)=MAT.mi_confusionB_BGwo;
            end
        end
    end
end
optWin=optWin(1:ii);
optWin_BGwo=optWin_BGwo(1:ii);

%% Plot the histograms of optimum window
figure(1)
subplot(2,1,1)
hist(optWin, 600)
%axis([0 610 0 300])
hold on
xlabel('WinSize of confusion matrix')
ylabel('nb cells')
title(sprintf('MiConfB nb of cells:%d', ii));
subplot(2,1,2)
hist(optWin_BGwo, 600)
%axis([0 610 0 300])
hold on
xlabel('WinSize of confusion matrix')
ylabel('nb cells')
title(sprintf('MiConfB without Backgound category, nb of cells: %d',ii))

WinSize=unique(optWin);
%% Plot the values of mi_conf for all cells depending on winsize
figure(2)
CO ={'b' ; 'g' ; 'r' ; 'c' ; 'm' ; 'y' ; 'k'}; 
for jj = 1:ii
    subplot(2,1,1)
    rn =randperm(7);
    plot(MI_confB(jj,:)/max(MI_confB(jj,:)), CO{rn(1)});
    set(gca,'XTickLabel',WinSize);
    axis([1 9 0 1.1])
    hold on
    subplot(2,1,2)
    plot(MI_confBBGwo(jj,:)/max(MI_confBBGwo(jj,:)), CO{rn(1)});
    axis([1 9 0 1.1])
    set(gca,'XTickLabel',WinSize);
    hold on
    %pause
end

hold off

figure(3)
plot(log2(MI_confB'))
hold on
set(gca,'XTickLabel',WinSize);
hold off