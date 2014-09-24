
cd /Users/elie/Documents/MATLAB/data/calmatfile
input_dir=pwd;
Subjects = dir(input_dir);
CallCat={ 'Ag' 'Be' 'DC' 'Di' 'LT' 'Ne' 'Te' 'Th' 'Wh' 'song' 'mlnoise'};
NCat=length(CallCat);

%creating the vector of all the type of comparison between vocalization hat
%exist
Ldsel=cumsum(1:(NCat));
Ldsel=Ldsel(end);
Compcat=cell(Ldsel,2);
Ncomp=0;
for vt1=1:(NCat)
    for vt2=vt1:NCat
        Ncomp=Ncomp+1;
        codevoc1=sum(double(char(CallCat(vt1))));
        codevoc2=sum(double(char(CallCat(vt2))));
        if codevoc1<=codevoc2
            Compcat{Ncomp,1}=[CallCat(vt1) CallCat(vt2)];
            Compcat{Ncomp,2}= codevoc1*codevoc2;
        elseif codevoc1>codevoc2
            Compcat{Ncomp,1}=[CallCat(vt2) CallCat(vt1)];
            Compcat{Ncomp,2}= codevoc1*codevoc2;
        end
    end
end


% Calculate the total number of units
Nunits=0;
for ss=1:length(Subjects)
    Subjects(ss).name=Indiv;
    if length(Indiv)==11
        Idir=fullfile(input_dir, Indiv);
        calfiles=dir(fullfile(Idir,'Cal*.mat'));
        LC=length(calfiles);
        Nunits = LC + Nunits;
    end
end

% Initialize matrices for average dprimes, StdDprimes, Dsel, AvgRcallcat,
% StdRcallcat
MDsel=zeros(Nunits,Ncomp);
MAvgDprime=zeros(Nunits, Ncomp);
MStdDprime=zeros(Nunits, Ncomp);
MAvgRcallcat=zeros(Nunits, NCat);
MStdRcallcat=zeros(Nunits, NCat);
UnitLRI=zeros(Nunits, NCat);
UnitSSI=zeros(Nunits, NCat);
UnitDprimeMeans=zeros(Nunits,1);
Unitnames=cell(Nunits);
UnitR2A=zeros(Nunits, 3);
UnitPValLRatio=zeros(Nunits, 2);


% Loop through units matfiles and fill in matrices
unit=0;
for ss=1:length(Subjects)
    Indiv= Subjects(ss).name;
    if length(Indiv)==11
        Idir=fullfile(input_dir, Indiv);
        calfiles=dir(fullfile(Idir,'Cal*nc.mat'));
        LC=length(calfiles);
        fprintf('Loading Matfiles of %s\n',Indiv);
        for cc=1:LC
            unit = unit + 1;
            calfile=fullfile(Idir, calfiles(cc).name);
            fprintf('Loading Matfile %s\n', calfiles(cc).name);
            Cal=load(calfile);
            
            % checking the type of vocalization category we have here
            outputI=1:NCat;
            inputI=zeros(1,length(Cal.StimType));
            ii=0;
            for vv=1:NCat
                ct = CallCat(vv);
                selector=find(strcmp(ct , Cal.StimType));
                if isempty(selector)
                    II=find(outputI==vv);
                    outputI(II)=[];
                elseif length(selector)==1
                    ii=ii+1;
                    inputI(ii)=selector;
                else
                    fprintf('could not correctly attribute call type');
                end
            end
            
            % filling matrices with results
            MAvgRcallcat(unit,outputI)= Cal.AvgRcallcat(inputI);
            MStdRcallcat(unit,outputI)= Cal.StdRcallcat(inputI);
            UnitLRI(unit,outputI)=Cal.LRIcat(inputI);
            UnitSSI(unit,outputI)=Cal.SSIcat(inputI);
            
            % checking the type of comparison between vocalization that was
            % made for Dsel calculations
            outputI_comp1=1:Ncomp;
            inputI_comp1=zeros(1,length(Cal.Dsel.comptype));
            ii=0;
            for vv=1:Ncomp
                ct = Compcat{vv,2};
                selector=find(cell2mat(Cal.Dsel.Comptype(:,2))==ct);
                if isempty(selector)
                    II=find(outputI_comp1==vv);
                    outputI_comp1(II)=[];
                elseif length(selector)==1
                    ii=ii+1;
                    inputI_comp1(ii)=selector;
                else
                    fprintf('could not correctly attribute comparison type');
                end
            end
            
            % filling matrix with results
            MDsel(unit,outputI_comp1) = cell2mat(Cal.Dsel.values(inputI_comp1));
            
            
            % checking the type of comparison between vocalization that was
            % made for AvgDprime calculations
            outputI_comp2=1:Ncomp;
            inputI_comp2=zeros(1,length(Cal.Comptype));
            ii=0;
            for vv=1:Ncomp
                ct = Compcat{vv,2};
                selector=find(cell2mat(Cal.Comptype(:,2))==ct);
                if isempty(selector)
                    II=find(outputI_comp2==vv);
                    outputI_comp2(II)=[];
                elseif length(selector)==1
                    ii=ii+1;
                    inputI_comp2(ii)=selector;
                else
                    fprintf('could not correctly attribute comparison type');
                end
            end
                
            
            MAvgDprime(unit,outputI_comp2)= cell2mat(Cal.AvgDprime(inputI_comp2));
            MStdDprime(unit,outputI_comp2)= cell2mat(Cal.StdDprime(inputI_comp2));
            
            UnitDprimeMeans(unit)=mean(MAvgDprime(unit,outputI_comp2));
            UnitR2A(unit, :)=Cal.R2A;
            UnitPValLRatio(unit, :)=Cal.PValLRatio(2:3);
            Unitnames{unit}=[Indiv calfiles(cc).name];
        end
    end
    fprintf('Done pulling data for %s\n',Subjects(ss).name);
end
%list of the different comparisons made
MComptype=Cal.Comptype(:,1);
             
 % Save matrices in a matfile
 if ismac()
        [status username] = system('who am i');
        if strcmp(strtok(username), 'frederictheunissen')
            if strncmp('/auto/fdata/solveig',stim_name, 19)
            elseif strncmp('/auto/fdata/julie',stim_name, 17)
                Dfilename = fullfile('/Users','frederictheunissen','Documents','Data','Julie','matfile',Res.subject,['Cal' Res.Site '.mat']);
            end
        elseif strcmp(strtok(username), 'elie')
            Dfilename = fullfile('/Users','elie','Documents','MATLAB','data','matfile','Dprimes', 'DprimeCalls.mat');
        end
else
    Dfilename=fullfile('/auto','k6','julie','matfile','Dprimes','DprimeCalls.mat');
end
save(Dfilename, 'M*', 'U*', 'Compcat', 'CallCat');
 
%% Play with the file ;-)
load('DprimeCalls.mat');
Nunits=length(Unitnames);

%% Set up the selector for voctype comparisons
Ag=zeros(Nunits,55);
Be = zeros(Nunits,55);
DC = zeros(Nunits,55);
Di = zeros(Nunits,55);
LT = zeros(Nunits,55);
Ne = zeros(Nunits,55);
Te = zeros(Nunits,55);
Th = zeros(Nunits,55);
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
%axis([0 2 0 150]);

%% Histogram of Dsel for all units
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
        subplot(2,LI/2,icomp)
        hist(MAvgDprimeS(:,IC));
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
 for vt=1:10
    subplot(2,5,vt)
    hist(MAvgRcallcat(:,vt));
    ytit=MComptype{(45+vt),1};
    title(char(strcat(ytit(1), ytit(2), '_', 'Mean spike rate /s')));
end