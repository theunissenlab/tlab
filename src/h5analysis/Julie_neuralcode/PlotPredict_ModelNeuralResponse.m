cd /Users/elie/Documents/MATLAB/data/matfile/CalDataBPCssmlnoise
input_dir=pwd;
Subjects = dir(input_dir);
Subj=3:7;
Legend= cell(3,1);
Legend{1}='Acoustic only';
Legend{2}='Semantic only';
Legend{3}='Acoustic + Semantic';
SignifModelCompAcSe= cell(5,1);
SignifModelCompAcTot= cell(5,1);
SignifModelCompSeTot=cell(5,1);
PValLRatioAcSe= cell(5,1);
PValLRatioAcTot=cell(5,1);
PValLRatioSeTot=cell(5,1);
R2AAc= cell(5,1);
R2ASe= cell(5,1);
R2ATot=cell(5,1);
for SS=1:5
    ss=Subj(SS);
    Idir=fullfile(input_dir, Subjects(ss).name);
    matfiles=dir(fullfile(Idir,'Cal*.mat'));
    LM=length(matfiles);
    SignifModelCompAcSe{SS}= zeros(LM,1);
    SignifModelCompAcTot{SS}= zeros(LM,1);
    SignifModelCompSeTot{SS}= zeros(LM,1);
    PValLRatioAcSe{SS}= zeros(LM,1);
    PValLRatioAcTot{SS}= zeros(LM,1);
    PValLRatioSeTot{SS}= zeros(LM,1);
    R2AAc{SS}= zeros(LM,1);
    R2ASe{SS}= zeros(LM,1);
    R2ATot{SS}= zeros(LM,1);
    fprintf('Reading Calfile of %s\n',Subjects(ss).name);
    for hh=1:LM
        MatfilePath=fullfile(Idir, matfiles(hh).name);
        fprintf('Reading %s\n', MatfilePath);
        Cal=load(MatfilePath);
        
        
        
        %reshape voctype I can supress this part if I run again
        %MatfilesConstruct_JEE (10/10/2012)
        nvoc=1;
DataSel=zeros(1,length(Cal.VocType));
voctype=Cal.VocType;
VT = zeros(1,9);
Ag = zeros(1,9);
Ag(1)=1;
Be = zeros(1,9);
Be(2)=1;
DC = zeros(1,9);
DC(3)=1;
Di = zeros(1,9);
Di(4)=1;
LT = zeros(1,9);
LT(5)=1;
Ne = zeros(1,9);
Ne(6)=1;
Te = zeros(1,9);
Te(7)=1;
Th = zeros(1,9);
Th(8)=1;
Song = zeros(1,9);
Song(9)=1;
for dd=1:length(voctype);
    if strcmp(voctype{dd}, 'Ag')
        VT = VT + Ag;
        DataSel(nvoc)=dd;
        nvoc=nvoc+1;
    elseif strcmp(voctype{dd}, 'Be')
        VT = VT + Be;
        DataSel(nvoc)=dd;
        nvoc=nvoc+1;
    elseif strcmp(voctype{dd}, 'DC')
        VT = VT + DC;
        DataSel(nvoc)=dd;
        nvoc=nvoc+1;
    elseif strcmp(voctype{dd}, 'Di')
        VT = VT + Di;
        DataSel(nvoc)=dd;
        nvoc=nvoc+1;
    elseif strcmp(voctype{dd}, 'LT')
        VT = VT + LT;
        DataSel(nvoc)=dd;
        nvoc=nvoc+1;       
    elseif strcmp(voctype{dd}, 'Ne')
        VT = VT + Ne;
        DataSel(nvoc)=dd;
        nvoc=nvoc+1;
    elseif strcmp(voctype{dd}, 'Te')
        VT = VT + Te;
        DataSel(nvoc)=dd;
        nvoc=nvoc+1;
    elseif strcmp(voctype{dd}, 'Th')
        VT = VT + Th;
        DataSel(nvoc)=dd;
        nvoc=nvoc+1;
    elseif strcmp(voctype{dd}, 'song')
        VT = VT + Song;
        DataSel(nvoc)=dd;
        nvoc=nvoc+1;
    end
end
VTMinNb=min(VT);
if (nvoc-1)==length(Cal.VocType)
    voc = Cal.VocType;
else
    DataSel=DataSel(DataSel>0);
    voc = Cal.VocType(DataSel);
end
        
        %Keep going with new list of voctype
        MAXY = max(cell2mat(Cal.ModelPredict));
        MAXX = max(Cal.NeuralResponse);
        SignifModelCompAcSe{SS}(hh)= Cal.SignifModelCompare(1);
        SignifModelCompAcTot{SS}(hh)= Cal.SignifModelCompare(2);
        SignifModelCompSeTot{SS}(hh)= Cal.SignifModelCompare(3);
        PValLRatioAcSe{SS}(hh)= Cal.PValLRatio(1);
        PValLRatioAcTot{SS}(hh)= Cal.PValLRatio(2);
        PValLRatioSeTot{SS}(hh)= Cal.PValLRatio(3);
        R2AAc{SS}(hh)= Cal.R2A(1);
        R2ASe{SS}(hh)= Cal.R2A(2);
        R2ATot{SS}(hh)= Cal.R2A(3);
        IndicesSelect=zeros(9,VTMinNb);
        LocalNeuralResponses =  zeros(9*VTMinNb,1);
        LocalModelPredict = cell(3,1);
        for mm = 1:3
            LocalModelPredict{mm} = zeros(9*VTMinNb,1);
        end
        LocalVoc = cell(9*VTMinNb,1);
        [C,IA,IC] = unique(voc);
        for ii = 1:9
            LocalDataSel=find(IC==ii);
            ID = randperm(length(LocalDataSel));
            Vocrand=LocalDataSel(ID);
            IndicesSelect(ii,1:VTMinNb) = Vocrand(1:VTMinNb);
            LocalNeuralResponses(((ii-1)*VTMinNb+1):ii*VTMinNb)=Cal.NeuralResponse(IndicesSelect(ii,:)); 
            for mm = 1:3
                LocalModelPredict{mm}(((ii-1)*VTMinNb+1):ii*VTMinNb) =  Cal.ModelPredict{mm}(IndicesSelect(ii,:));
            end
            LocalVoc(((ii-1)*VTMinNb+1):ii*VTMinNb) = voc(IndicesSelect(ii,:));
        end
       
        for jj=1:3
            figure(1)
            subplot(1,3,jj);
            gscatter(LocalNeuralResponses, LocalModelPredict{jj}, LocalVoc, 'mgcbrkyyr', '......d.d',[20 20 20 20 20 20 10 20 10]);
            %gscatter(LocalNeuralResponses, LocalModelPredict{jj}, LocalVoc, 'rkmgcbryy', '......dd.',[20 20 20 20 20 20 10 10 20]);
            title(sprintf('%s Adjusted R^2 = %.2f', Legend{jj}, Cal.R2A(jj,1)));
            MAX=max(MAXX,MAXY);
            axis([0 MAX 0 MAX]);
            hold on
            plot(0:MAX/10:MAX,0:MAX/10:MAX, 'k');
            hold off
        end
       % subplot(1,4,4);
        %plot(1:3, Cal.R2A);
        %axis([1 3 0 1]);
        fprintf('**Acoustic only: %f\nSemantic only: %f\nAcoustic + Semantic: %f\n', Cal.R2A(1,1), Cal.R2A(2,1), Cal.R2A(3,1));
        fprintf('**Acoustic vs Semantic: %f\n Acoustic vs Tot: %f\nSemantic vs Tot: %f\n', Cal.PValLRatio(1), Cal.PValLRatio(2), Cal.PValLRatio(3));
        pause
    end
end


    figure(2)
    subplot(3,1,1);
    hist(cell2mat(R2AAc),20);
    mean(cell2mat(R2AAc))
    std(cell2mat(R2AAc))
    axis([0 1 0 25])
    subplot(3,1,2);
    hist(cell2mat(R2ASe),20);
    mean(cell2mat(R2ASe))
    std(cell2mat(R2ASe))
    axis([0 1 0 25])
    subplot(3,1,3);
    hist(cell2mat(R2ATot),20);
    mean(cell2mat(R2ATot))
    std(cell2mat(R2ATot))
    axis([0 1 0 25])
    
    figure(5)
    subplot(3,1,1);
    R2ratioSeAc=cell2mat(R2ASe)./cell2mat(R2AAc);
    sum(R2ratioSeAc>1)
    hist(R2ratioSeAc);
    %axis([-2.5 2.5 0 100])
    axis([0 4.5 0 140])
    subplot(3,1,2);
    R2ratioTotSe=cell2mat(R2ATot)./cell2mat(R2ASe); 
    hist(R2ratioTotSe);
    %axis([-2.5 2.5 0 100])
    axis([0 4.5 0 140])
    subplot(3,1,3);
    R2ratioTotAc=cell2mat(R2ATot)./cell2mat(R2AAc); 
    hist(R2ratioTotAc);
    %axis([-2.5 2.5 0 100])
     axis([0 4.5 0 140])
    
    
    figure(3)
    subplot(1,3,1);
    hist(cell2mat(SignifModelCompAcSe));
    subplot(1,3,2);
    hist(cell2mat(SignifModelCompAcTot));
    subplot(1,3,3);
    hist(cell2mat(SignifModelCompSeTot));
    
    figure(4)
    subplot(1,3,1);
    hist(cell2mat(PValLRatioAcSe));
    axis([0 1 0 180])
    subplot(1,3,2);
    hist(cell2mat(PValLRatioAcTot));
    %axis([0 1 0 180])
    subplot(1,3,3);
    hist(cell2mat(PValLRatioSeTot));
    %axis([0 1 0 180])