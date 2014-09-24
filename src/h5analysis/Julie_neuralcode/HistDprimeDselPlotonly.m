%% Play with the file ;-)
load('AllUnits.mat');
Nunits=length(Unitnames);
CallCatI=[1:8 10];
CallCat=CallCat(CallCatI); %get rid of Wh and mln for this analysis because not enough cells for Whi (5 only)

NStimType=length(CallCat);

%% Calculate good LRI
UnitLRI2=NaN(Nunits,NStimType);
for uu=1:Nunits
    Othermean=zeros(NStimType,1);
    Selfmean=zeros(NStimType,1);
    for st= 1:NStimType
        cc=CallCatI(st);
        Othermean_sel=MAvgRcallcat(uu,CallCatI);
        Othermean_sel(st)=[]; %supress the category we are looking at from others
        Othermean(st)=nanmean(Othermean_sel);
        Selfmean(st)=MAvgRcallcat(uu,cc);
    end
    Selfmean(find(Selfmean==0))=min(Selfmean(find(isfinite(log2(Selfmean)))))/100; % supress zero values from that vector
    UnitLRI2(uu,:)=log2(Selfmean./Othermean);
end

%% Set up the selector for voctype comparisons % NOT USEFULL ANYMORE?
Ag=zeros(Nunits,55);
Be = zeros(Nunits,55);
DC = zeros(Nunits,55);
Di = zeros(Nunits,55);
LT = zeros(Nunits,55);
Ne = zeros(Nunits,55);
Te = zeros(Nunits,55);
Th = zeros(Nunits,55);
Wh = zeros(Nunits, 55);
mlnoise = zeros(Nunits,55);
song = zeros(Nunits,55);
Sign=[1 -1]; % if voctype i is in second column then the result is given for j-i so results should be * by -1 to get the results for i-j
for ll=1:length(MComptype)
    N=char(MComptype{ll});
        for ii=1:2
            if strncmp(N(ii,:),'Ag',2)
                Ag(:,ll)=Sign(ii);
            elseif strncmp(N(ii,:),'Be',2)
                Be(:,ll)=Sign(ii);
            elseif strncmp(N(ii,:),'DC',2)
                DC(:,ll)=Sign(ii);
            elseif strncmp(N(ii,:),'Di',2)
                Di(:,ll)=Sign(ii);
            elseif strncmp(N(ii,:),'LT',2)
                LT(:,ll)=Sign(ii);
            elseif strncmp(N(ii,:),'Ne',2)
                Ne(:,ll)=Sign(ii);
            elseif strncmp(N(ii,:),'Te',2)  
                Te(:,ll)=Sign(ii);
            elseif strncmp(N(ii,:),'Th',2)
                Th(:,ll)=Sign(ii);
            elseif strncmp(N(ii,:),'mlnoise',2)
                mlnoise(:,ll)=Sign(ii);
            elseif strncmp(N(ii,:),'song',2)
                song(:,ll)=Sign(ii);
            end
        end
end
VocSelector={Ag, Be, DC, Di, LT, Ne, Te, Th, mlnoise, song};


%% Plot mean Dprime over all stim category of all units
figure(1)
hist(UnitDprimeMeans,100); %to track cells that does not respond differently to any sound
title('Mean Dprime over all call categories');
vline(0.3)
sum(UnitDprimeMeans<0.2)
%axis([0 2 0 150]);

%% Histogram of Dsel for all units This migt not be working anymore because I sligtly changed data format
%Be aware of squewed distribution whenever DC or beggings are involved
% Dsel is calculated as (mean(spike rate stimulus category 1)-mean(spike
% rate stimulus category 2))/(mean(spike rate stimulus category 1)+mean(spike
% rate stimulus category 2)). the legend respect the order in which Dsel
% was calculated for each hisogram eg: AgBe -> Ag-Be
figure(2)
for icomp=1:45
    subplot(5,9,icomp)
    hist(MDsel(:,icomp));
    ytit=MComptype{icomp};
    title(char(strcat(ytit(1), ytit(2))));
end

%% Histograms of Dsel for all units for all comparisons per voc category
%% This migt not be working anymore because I sligtly changed data format
%same thing as above but sort by stim category
% for each graph, the value are given by (stim being compared)-(all the other stims)
% for instance Ag-Be, Ag-DC, Ag-Di... for the first graph even if the
% legend is not written in the good way
figure(3)
for vt=1:10
    Sel=VocSelector{vt};
    MDselS=MDsel.*Sel(:,1:45);
    Indsel=find(Sel(1,1:45));
    LI=length(Indsel);
    for icomp=1:LI
        IC=Indsel(icomp);
        subplot(2,ceil(LI/2),icomp)
        hist(MDselS(:,IC));
        ytit=MComptype{IC};
        title(char(strcat(ytit(1), ytit(2))));
    end
    pause
end



%% Histogram of Average Dprime for all units
% the legend respect the order in which Dsel
% was calculated for each hisogram eg: AgBe -> Ag-Be
figure(4)
for icomp=1:55
    subplot(5,11,icomp)
    hist(MAvgDprime(:,icomp));
    %axis([5 20 0 25]);
    ytit=MComptype{icomp};
    title(char(strcat(ytit(1), ytit(2), '_', 'Dprime')));
end

%% Histograms of average Dprime for all units per voc category
%same thing as above but sort by stim category
% for each graph, the value are given by (stim being compared)-(all the other stims)
% for instance Ag-Be, Ag-DC, Ag-Di... for the first graph even if the
% legend is not written in the good way
 figure(5)
 for vt=1:10
    Sel=VocSelector{vt};
    MAvgDprimeS=MAvgDprime.*Sel;
    Indsel=find(Sel(1,:));
    LI=length(Indsel);
    for icomp=1:LI
        IC=Indsel(icomp);
        subplot(1,LI,icomp)
        hist(MAvgDprimeS(:,IC));
        axis([-10 10 0 50])
        ytit=MComptype{IC};
        title(char(strcat(ytit(1), ytit(2),'  ', 'Mean spike rate /s')));
    end
    pause
 end
 
 %% Barplot per unit of Dprime value per voc category
% for each graph, the value are given by (stim being compared)-(all the other stims)
% for instance Ag-Be, Ag-DC, Ag-Di... for the first graph even if the
% legend is not written in the good way
 figure(6)
 for unit=1:Nunits
     for vt=1:10
        Sel=VocSelector{vt};
        MAvgDprimeS=MAvgDprime.*Sel;
        MStdDprimeS=MStdDprime.*Sel;
        Indsel=find(Sel(1,:));
        LI=length(Indsel);
        bar(MAvgDprimeS(unit,Indsel));
        Xtit='';
        for xx=1:LI
            ind=Indsel(xx);
            xtit=MComptype{ind};
            Xtit=char(strcat(Xtit, xtit(1), xtit(2),'--'));
        end
        title(char(strcat('Mean spike rate /s', char(Unitnames{unit}))));
        xlabel(Xtit);
     pause
     end
 end
 
 %% Histogram of average spike rate for all units
 figure(7)
 for vt=1:NStimType
     cc=CallCatI(vt);
    subplot(3,4,vt)
    hist(MAvgRcallcat(:,cc),100);
    axis([0 600 0 150]);
    Moy=nanmean(MAvgRcallcat(:,cc));
    SD=nanstd(MAvgRcallcat(:,cc));
    vline(Moy);
    vline(Moy+SD, 'g');
    vline(Moy-SD, 'g');
    title(char(strcat(CallCat(vt), '_', 'Mean spike rate /s')));
    text(200, 100, sprintf('Mean =%0.2f +/- %0.2f', Moy, SD));
 end
 
 %% mean LRI  and mean SSI for each call category
 figure(8)
 subplot(2,1,1)
 Mean_LRI=NaN(NStimType,1);
 Std_LRI=NaN(NStimType,1);
 for cc=1:NStimType
     LRI_finite=find(isfinite(UnitLRI2(:,cc)));
     Mean_LRI(cc)=nanmean(abs(UnitLRI2(LRI_finite,cc)));
     Std_LRI(cc)=nanstd(abs(UnitLRI2(LRI_finite,cc)));
 end
 bar(Mean_LRI);
 hold on
 errorbar(Mean_LRI, Std_LRI, 'b.');
 set(gca(), 'Xtick', 1:NStimType);
 set(gca(), 'XTickLabel', CallCat(1:NStimType));
 
 subplot(2,1,2)
 Mean_SSI=NaN(NStimType,1);
 Std_SSI=NaN(NStimType,1);
 for st=1:NStimType
     cc=CallCatI(st);
     SSI_finite=find(isfinite(UnitSSI(:,cc)));
     Mean_SSI(cc)=nanmean(abs(UnitSSI(SSI_finite,cc)));
     Std_SSI(cc)=nanstd(abs(UnitSSI(SSI_finite,cc)));
 end
 bar(Mean_SSI);
 hold on
 errorbar(Mean_SSI, Std_SSI, 'b.');
 set(gca(), 'Xtick', 1:NStimType);
 set(gca(), 'XTickLabel', CallCat(1:NStimType));
 hline(1/3)
 %% representing, for each call category the percentage of multiunits that have a LRI (Log Ratio Index) above a certain threshold
 
 figure(9)
 sp =0;
 Thresholds=[-1 -2 -3 1 2 3];
 for Thr=1:6
     thr=Thresholds(Thr);
     PercCell_local=NaN(NStimType,1);
     for cc=1:NStimType
        UnitsIdx=~isnan(UnitLRI2(:,cc));
        NbUnits_local=sum(UnitsIdx);
        UnitsIdx=find(UnitsIdx);
        UnitLRI_local=UnitLRI2(UnitsIdx,cc);
        if thr>0
            PercCell_local(cc)=100*sum(UnitLRI_local>=thr)/NbUnits_local;
        elseif thr<0
            PercCell_local(cc)=100*sum(UnitLRI_local<=thr)/NbUnits_local;
        end
     end
     sp=sp+1;
     subplot(2,3,sp)
     bar(PercCell_local);
     if thr>0
        title(sprintf('Percentage of multiunits with LRI>=%d', thr));
     elseif thr<0
         title(sprintf('Percentage of multiunits with LRI<=%d', thr));
     end
     set(gca(), 'Xtick', 1:NStimType);
     set(gca(), 'XTickLabel', CallCat(1:NStimType));
 end
 
 %% representing the distribution of multiunits depending of the number of LRI be 1.
 %Using SSI with a threshold of 1/3 to find cell that discriminate at least
 %one category
 figure(10)
 nbHighLRI=zeros(Nunits,1);
 nbHighSSI=zeros(Nunits,1);
 for uu=1:Nunits
    nbHighLRI(uu)=sum(abs(UnitLRI2(uu,:))>=1);
    nbHighSSI(uu)=sum(abs(UnitSSI(uu,:))>=1/3);
 end
 
 hist(nbHighLRI, max(nbHighLRI)+1)
 axis([-0.5 max(nbHighLRI)+1 0 Nunits])
 H=hist(nbHighLRI, max(nbHighLRI)+1)
 LRIthr=0:max(nbHighLRI);
 for tt=1:length(H)
     fprintf('%f%% of units present %d abs(LRI) values above or equal to 1\n', 100*H(tt)/Nunits, LRIthr(tt));
 end
 HighLRIUnits=find(nbHighLRI>0); %indices of the units that have at least 1 category of vocalizations with LRI>1 or LRI<-1
 LowLRIUnits=find(nbHighLRI==0); %indices of units that are not considered as discriminative of vocalization categories according to threshold
 HighSSIUnits=find(nbHighSSI>0); %should be the same as HighLRIUnits
 LowSSIUnits=find(nbHighSSI==0); %should be the same as LowLRIUnits
 
 
 %% represent R2A values of acoustic model
 
 Models={'Acoustic' 'Semantic' 'Acoustic+Semantic'};
 for rr=1:3
    figure(11)
    subplot(2,3,rr)
    hist(UnitR2A(HighLRIUnits, rr))
    vline(mean(UnitR2A(HighLRIUnits, rr)));
    axis([-0.5 1 0 40])
    title(strcat(Models{rr},' Units with nb(LRI>1)>1'));
    text(-0.5, 35, sprintf('mean R2A %f',mean(UnitR2A(HighLRIUnits, rr)))); 
    subplot(2,3,rr+3)
    hist(UnitR2A(LowLRIUnits,rr))
    vline(mean(UnitR2A(LowLRIUnits,rr)));
    axis([-0.5 1 0 200])
    title(strcat(Models{rr}, ' Units with all LRI<1'));
    text(0.1, 175,sprintf('mean R2A %f',mean(UnitR2A(LowLRIUnits, rr))));
    
    figure(12)
    subplot(1,4,rr)
    plot(nbHighLRI, UnitR2A(:,rr),'bx')
    title(strcat(Models{rr},' R2A against nb(LRI>1)'))
    axis([0 8 -0.1 1])
 end
 subplot(1,4,4)
 plot(nbHighLRI, UnitR2A(:,3)-UnitR2A(:,1), 'bx');
 axis([0 8 -0.1 1])
   
 %% plot std LRI for each multiunit wit R2A
 STDLRI=zeros(Nunits,1);
 for uu=1:Nunits
     FiniteVal=find(isfinite(UnitLRI2(uu,:)));
    STDLRI(uu)=nanstd(UnitLRI2(uu,FiniteVal));
 end
 
 figure(13)
 plot(nbHighLRI, STDLRI, 'bx'); %good correlation between the two metrics
 axis([0 10 0 3])
 xlabel('nb LRI above 1')
 ylabel('std LRI')
 mdl=LinearModel.fit(nbHighLRI,STDLRI);
 R2A=mdl.Rsquared.Adjusted;
 [R2, P2]=corrcoef([nbHighLRI STDLRI])
 text(0.5, 2, sprintf('R2a=%f', R2A));
 Intercept=mdl.Coefficients.Estimate(1);
 slope=mdl.Coefficients.Estimate(2);
 line('XData', [0 10], 'YData', [Intercept 10*slope+Intercept], 'LineStyle', '-', ...
    'LineWidth', 2, 'Color','m');
 
 figure(14)
 for rr=1:3
    subplot(1,5,rr)
    plot(STDLRI, UnitR2A(:,rr),'bx')
    title(strcat(Models{rr},' R2A against std(LRI)'))
    xlabel('std(LRI)')
    ylabel(strcat('R2a ',Models{rr},' model'))
    mdl=LinearModel.fit(STDLRI, UnitR2A(:,rr));
    R2A=mdl.Rsquared.Adjusted;
    [R2, P2]=corrcoef([STDLRI UnitR2A(:,rr)])
    text(0.1, 0.95, sprintf('R2a = %f', R2A));
    axis([0 2 -0.1 1])
    Intercept=mdl.Coefficients.Estimate(1);
    slope=mdl.Coefficients.Estimate(2);
    line('XData', [0 2], 'YData', [Intercept 2*slope+Intercept], 'LineStyle', '-', ...
    'LineWidth', 2, 'Color','m');
 end
 subplot(1,5,4)
 plot(STDLRI, UnitR2A(:,3)-UnitR2A(:,1), 'bx');
 title('R2A(A+S)-R2A(A) against std(LRI)')
 xlabel('std(LRI)')
 ylabel('R2a(Acoustic + Semantic) - R2a(Acoustic)');
 DiffR2A=UnitR2A(:,3)-UnitR2A(:,1);
 mdl2=LinearModel.fit(STDLRI, DiffR2A);
 R2A2=mdl2.Rsquared.Adjusted;
 text(0.1, 0.19, sprintf('R2a=%f', R2A2));
 axis([0 2 -0.1 0.2])
 Intercept=mdl2.Coefficients.Estimate(1);
 slope=mdl2.Coefficients.Estimate(2);
 line('XData', [0 2], 'YData', [Intercept 2*slope+Intercept], 'LineStyle', '-', ...
    'LineWidth', 2, 'Color','m');
 

 subplot(1,5,5)
 SignifI=find((UnitPValLRatio(:,1)<0.05));
 plot(STDLRI(SignifI), DiffR2A(SignifI), 'bx');
 mdl2=LinearModel.fit(STDLRI(SignifI), DiffR2A(SignifI));
 R2A2=mdl2.Rsquared.Adjusted;
 title('signif R2A(A+S)-R2A(A) against std(LRI)')
 xlabel('std(LRI)')
 ylabel('signif R2a(Acoustic + Semantic) - R2a(Acoustic)')
 axis([0 2 -0.1 0.2])
 text(0.1, 0.19, sprintf('R2a=%f', R2A2));
 Intercept=mdl2.Coefficients.Estimate(1);
 slope=mdl2.Coefficients.Estimate(2);
 line('XData', [0 2], 'YData', [Intercept 2*slope+Intercept], 'LineStyle', '-', ...
    'LineWidth', 2, 'Color','m');

%% Exploring the discriminative cells, those with at least one LRI>1 or one SSI>1/3
%% What category of calls is the most discriminated by neurons

figure(15)
 sp =0;
 Thresholds=[-1 -2 -3 1 2 3];
 %Thresholds=[-7/9 -3/5 -1/3 1/3 3/5 7/9];
 for Thr=1:6
     thr=Thresholds(Thr);
     PercCell_local=NaN(NStimType,1);
     for cc=1:NStimType
        UnitLRI_local=UnitLRI2(HighLRIUnits,cc);
        UnitLRI_local=UnitLRI_local(find(~isnan(UnitLRI_local)));
        NbUnits_local=length(UnitLRI_local);
        if thr>0
            PercCell_local(cc)=100*sum(UnitLRI_local>=thr)/NbUnits_local;
        elseif thr<0
            PercCell_local(cc)=100*sum(UnitLRI_local<=thr)/NbUnits_local;
        end
     end
     sp=sp+1;
     subplot(2,3,sp)
     bar(PercCell_local);
     if thr>0
        title(sprintf('Percentage of multiunits with LRI>=%d', thr));
     elseif thr<0
         title(sprintf('Percentage of multiunits with LRI<=%d', thr));
     end
     set(gca(), 'Xtick', 1:NStimType);
     set(gca(), 'XTickLabel', CallCat(1:NStimType));
 end
 PercCell=NaN(NStimType,3);
 for cc=1:NStimType
    UnitLRI_local=UnitLRI2(HighLRIUnits,cc);
    UnitLRI_local=UnitLRI_local(find(~isnan(UnitLRI_local)));
    NbUnits_local=length(UnitLRI_local);
    PercCell(cc,1)=100*sum(UnitLRI_local<=-1)/NbUnits_local;
    PercCell(cc,2)=100*sum(abs(UnitLRI_local)<1)/NbUnits_local;
    PercCell(cc,3)=100*sum(UnitLRI_local>=1)/NbUnits_local;
 end
 
 figure(16)
 sp =0;
 Thresholds=[1 2 3 4 5 6 7];
 %Thresholds=[-7/9 -3/5 -1/3 1/3 3/5 7/9];
 for thr=1:7
     PercCell_local=NaN(NStimType,2);
     for cc=1:NStimType
       
        UnitLRI_local=UnitLRI2(nbHighLRI==thr,cc);
        UnitLRI_local=UnitLRI_local(find(~isnan(UnitLRI_local)));
        NbUnits_local=length(UnitLRI_local);
        PercCell_local(cc,1)=100*sum(UnitLRI_local>=1)/NbUnits_local;
        PercCell_local(cc,2)=100*sum(UnitLRI_local<=-1)/NbUnits_local;
        
        UnitLRI_cc=UnitLRI2(:,cc); 
        NbUNITS=sum(abs(UnitLRI_cc)>0);
        HighUnit_cc=find(abs(UnitLRI_cc)>0);
        PerCell(thr,cc)=100*sum(nbHighLRI(HighUnit_cc)==thr)/NbUNITS;
     end
     sp=sp+1;
     subplot(7,2,sp)
     bar(PercCell_local(:,1));
     title(sprintf('Percentage of multiunits with nb(LRI>=1)=%d', thr));
     set(gca(), 'Xtick', 1:NStimType);
     set(gca(), 'XTickLabel', CallCat(1:NStimType));
     sp=sp+1;
     subplot(7,2,sp)
     bar(PercCell_local(:,2));
     title(sprintf('Percentage of multiunits with nb(LRI<=-1)=%d', thr));
     set(gca(), 'Xtick', 1:NStimType);
     set(gca(), 'XTickLabel', CallCat(1:NStimType));
 end
 
 PercCell=NaN(7, NStimType);
 for cc=1:NStimType
     UnitLRI_cc=UnitLRI2(:,cc); 
     NbUNITS=sum(abs(UnitLRI_cc)>=1);
     HighUnit_cc=find(abs(UnitLRI_cc)>=1);
     nbHighLRI_cc=nbHighLRI(HighUnit_cc);
     for thr=1:7
         PercCell(thr,cc)=100*sum(nbHighLRI_cc==thr)/NbUNITS;
     end
 end
     
       
 
 figure(17)
 hist(nbHighLRI(find(nbHighLRI>=1)), max(nbHighLRI))
 H=hist(nbHighLRI(find(nbHighLRI>=1)), max(nbHighLRI))
 LRIthr=1:max(nbHighLRI);
 for tt=1:length(H)
     fprintf('%f%% of units present %d abs(LRI) values above or equal to 1\n', 100*H(tt)/sum(nbHighLRI>=1), LRIthr(tt));
 end
 
 %% plot std AvgRCallCat for HighLRI or LowLRI multiunit wit R2A
 HighLRIUnits=find(nbHighLRI>0); %indices of the units that have at least 1 category of vocalizations with LRI>1 or LRI<-1
 LowLRIUnits=find(nbHighLRI==0);
 UnitCategories={HighLRIUnits, LowLRIUnits};
 for catu=1:2;
    UnitCat=UnitCategories{catu};
    local_Nunits=length(UnitCat);
    STDSR=zeros(local_Nunits,1);
    for ii=1:local_Nunits
        uu=UnitCat(ii);
        STDSR(ii)=nanstd(MAvgRcallcat(uu,1:9));
    end
 
    figure(17+catu)
    plot(nbHighLRI(UnitCat), STDSR, 'bx'); %no correlation between the two metrics
    %axis([0 10 0 3])
    xlabel('nb LRI above 1')
    ylabel('std mean Spike rate per category')
    mdl=LinearModel.fit(nbHighLRI(UnitCat),STDSR);
    R2A=mdl.Rsquared.Adjusted;
    [R2, P2]=corrcoef([nbHighLRI(UnitCat) STDSR])
    text(1.5, 59, sprintf('R2a=%f', R2A));
    Intercept=mdl.Coefficients.Estimate(1);
    slope=mdl.Coefficients.Estimate(2);
 line('XData', [0 10], 'YData', [Intercept 10*slope+Intercept], 'LineStyle', '-', ...
    'LineWidth', 2, 'Color','m');
 
    figure(19+catu)
    for rr=1:3
        subplot(1,5,rr)
        plot(STDSR, UnitR2A(UnitCat,rr),'bx')
        title(strcat(Models{rr},' R2A against std(Spike rate)'))
        xlabel('std(SR)')
        ylabel(strcat('R2a ',Models{rr},' model'))
        mdl=LinearModel.fit(STDSR, UnitR2A(UnitCat,rr));
        R2A=mdl.Rsquared.Adjusted;
        [R2, P2]=corrcoef([STDSR UnitR2A(UnitCat,rr)])
        text(1, 0.9, sprintf('R2a = %f', R2A));
        axis([0 70 0 1])
        Intercept=mdl.Coefficients.Estimate(1);
        slope=mdl.Coefficients.Estimate(2);
        line('XData', [0 70], 'YData', [Intercept 70*slope+Intercept], 'LineStyle', '-', ...
    'LineWidth', 2, 'Color','m');
    end
    subplot(1,5,4)
    plot(STDSR, UnitR2A(UnitCat,3)-UnitR2A(UnitCat,1), 'bx');
    title('R2A(A+S)-R2A(A) against std(SR)')
    xlabel('std(SR)')
    ylabel('R2a(Acoustic + Semantic) - R2a(Acoustic)');
    DiffR2A=UnitR2A(UnitCat,3)-UnitR2A(UnitCat,1);
    mdl2=LinearModel.fit(STDSR, DiffR2A);
    R2A2=mdl2.Rsquared.Adjusted;
    text(1, 0.19, sprintf('R2a=%f', R2A2));
    axis([0 70 0 0.2])
    Intercept=mdl2.Coefficients.Estimate(1);
    slope=mdl2.Coefficients.Estimate(2);
    line('XData', [0 70], 'YData', [Intercept 70*slope+Intercept], 'LineStyle', '-', ...
    'LineWidth', 2, 'Color','m');
 

    subplot(1,5,5)
    SignifI=find((UnitPValLRatio(UnitCat,1)<0.05));
    plot(STDSR(SignifI), DiffR2A(SignifI), 'bx');
    mdl2=LinearModel.fit(STDSR(SignifI), DiffR2A(SignifI));
    R2A2=mdl2.Rsquared.Adjusted;
    title('signif R2A(A+S)-R2A(A) against std(SR)')
    xlabel('std(SR)')
    ylabel('signif R2a(Acoustic + Semantic) - R2a(Acoustic)')
    axis([0 70 0 0.2])
    text(1, 0.19, sprintf('R2a=%f', R2A2));
    Intercept=mdl2.Coefficients.Estimate(1);
    slope=mdl2.Coefficients.Estimate(2);
    line('XData', [0 70], 'YData', [Intercept 70*slope+Intercept], 'LineStyle', '-', ...
    'LineWidth', 2, 'Color','m');
 end

 %% percentage of units that are well explained by the acoustic model
 Perc=sum(UnitPvalue<0.01)/Nunits*100;
 for mm=1:3
     MeanMM(mm)=mean(UnitR2A(find(UnitPvalue(:,mm)<0.01),mm));
     STDMM(mm)=std(UnitR2A(find(UnitPvalue(:,mm)<0.01),mm));
     fprintf('%0.2f%% units are significantly well explained by the acoustic model with R2a=%0.2f+/-%0.2f\n',Perc(mm),MeanMM(mm),STDMM(mm));
 end
 for mm=1:3
     MeanMM(mm)=mean(UnitR2A(:,mm));
     STDMM(mm)=std(UnitR2A(:,mm));
     fprintf('models present R2a=%0.2f+/-%0.2f\n',MeanMM(mm),STDMM(mm));
 end
 NunitsHigh=length(HighLRIUnits)
 NunitsLow=length(LowLRIUnits)
 UnitR2AHigh=UnitR2A(HighLRIUnits,:); %indices of the units that have at least 1 category of vocalizations with LRI>1 or LRI<-1
 UnitR2ALow=UnitR2A(LowLRIUnits,:);
 UnitPvalueHigh=UnitPvalue(HighLRIUnits,:);
 UnitPvalueLow=UnitPvalue(LowLRIUnits,:);
 PercHigh=sum(UnitPvalueHigh<0.01)/NunitsHigh*100;
 PercLow=sum(UnitPvalueLow<0.01)/NunitsLow*100;
 
 for mm=1:3
     MeanMM(mm)=mean(UnitR2AHigh(find(UnitPvalueHigh(:,mm)<0.01),mm));
     STDMM(mm)=std(UnitR2AHigh(find(UnitPvalueHigh(:,mm)<0.01),mm));
     fprintf('%0.2f%% high units are significantly well explained by the acoustic model with R2a=%0.2f+/-%0.2f\n',PercHigh(mm),MeanMM(mm),STDMM(mm));
 end
 
 for mm=1:3
     MeanMM(mm)=mean(UnitR2ALow(find(UnitPvalueLow(:,mm)<0.01),mm));
     STDMM(mm)=std(UnitR2ALow(find(UnitPvalueLow(:,mm)<0.01),mm));
     fprintf('%0.2f%% low units are significantly well explained by the acoustic model with R2a=%0.2f+/-%0.2f\n',PercLow(mm),MeanMM(mm),STDMM(mm));
 end
 
 %% Find a cell that respond strongly to DC only
 Ind1=find(nbHighLRI==1);
 Ind2=find(UnitLRI2(:,3)>=1);
 Indices=intersect(Ind1, Ind2);
 max(UnitLRI2(Indices,3));
 MAX=max(UnitLRI2(Indices,3))
 UNDC=find(UnitLRI2(:,3)==MAX)
 Uname=cell2mat(Unitnames(UNDC))
 h5PathDC=strcat('/Users/elie/Documents/MATLAB/data/h5/', Uname(1:11), '/', Uname(15:end-7), '.h5');
 inputargDC = {'number', '115'};
 
 inputargSo = {'number', '147'}
 inputargBE = {'number', '56-58'};
 h5_plot_select_spikes(h5PathDC, 'Call1', 1, inputargDC)
 h5_plot_select_spikes(h5PathDC, 'Call1', 1, inputargSo)
 figure(3)
    H=bar(MAvgRcallcat(UNDC,CallCatI));
    hold on
    h2=errorbar(MAvgRcallcat(UNDC,CallCatI), MStdRcallcat(UNDC,CallCatI), 'b.');
    set(gca(), 'Xtick', 1:length(CallCat));
    set(gca(), 'XTickLabel', {'Ag' 'Qu' 'DC' 'De' 'LT' 'Ni' 'Te' 'Th' 'Ch'});
 
 %% Find a cell that don't respond to Be only
 Ind1=find(nbHighLRI==1);
 Ind2=find(UnitLRI2(:,2)<=-1);
 Indices=intersect(Ind1, Ind2);
 min(UnitLRI2(Indices,2));
 MAX=min(UnitLRI2(Indices,2))
 UNBe=find(UnitLRI2(:,2)==MAX)
 Indices(find(Indices==UNBe))=[];
 MAX=min(UnitLRI2(Indices,2))
 UNBe=find(UnitLRI2(:,2)==MAX)
 Indices(find(Indices==UNBe))=[];
 MAX=min(UnitLRI2(Indices,2))
 UNBe=find(UnitLRI2(:,2)==MAX)
 
 %UNBe=146;
 Uname=cell2mat(Unitnames(UNBe))
 h5PathBe=strcat('/Users/elie/Documents/MATLAB/data/h5/', Uname(1:11), '/', Uname(15:end-7), '.h5');
 inputargDC = {'number', '115'};
 inputargSo = {'number', '78'}
 inputargBE = {'number', '62'};
 h5_plot_select_spikes(h5PathBe, 'Call1', 1, inputargBE)
 h5_plot_select_spikes(h5PathBe, 'Call1', 1, inputargSo)
  %h5_plot_select_spikes(h5PathBe, 'Call1', 1, inputargDC)
  
  figure(6)
    H=bar(MAvgRcallcat(UNBe,CallCatI));
    hold on
    h2=errorbar(MAvgRcallcat(UNBe,CallCatI), MStdRcallcat(UNBe,CallCatI), 'b.');
    set(gca(), 'Xtick', 1:length(CallCat));
    set(gca(), 'XTickLabel', {'Ag' 'Qu' 'DC' 'De' 'LT' 'Ni' 'Te' 'Th' 'Ch'});
   
    %%  Non selective cell
    Uname=cell2mat(Unitnames(136))
    h5PathNS=strcat('/Users/elie/Documents/MATLAB/data/h5/', Uname(1:11), '/', Uname(15:end-7), '.h5');
 inputargDC = {'number', '98'};
 inputargSo = {'number', '78'}
 inputargBE = {'number', '62'};
 h5_plot_select_spikes(h5PathNS, 'Call2', 1, inputargDC)
 h5_plot_select_spikes(h5PathNS, 'Call2', 1, inputargSo)
 figure(9)
    H=bar(MAvgRcallcat(136,CallCatI));
    hold on
    h2=errorbar(MAvgRcallcat(136,CallCatI), MStdRcallcat(136,CallCatI), 'b.');
    set(gca(), 'Xtick', 1:length(CallCat));
    set(gca(), 'XTickLabel', {'Ag' 'Qu' 'DC' 'De' 'LT' 'Ni' 'Te' 'Th' 'Ch'});
 
 