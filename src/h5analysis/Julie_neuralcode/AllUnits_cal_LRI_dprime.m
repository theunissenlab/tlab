function AllUnits_cal_LRI_dprime(UT)
%note that most of the analysis does not take whines into account ...
...(only average spike rate does take it into account) because
%we have data for only one bird for now. This should be changed if we have
%more whine datas

%cd /Users/elie/Documents/MATLAB/data/matfile
cd /auto/k6/julie/matfile
CallCatRate={ 'Ag' 'Be' 'DC' 'Di' 'LT' 'Ne' 'Te' 'Th' 'song' 'Wh' 'mlnoise'};
CallCatDP={'Ag' 'Be' 'DC' 'Di' 'LT' 'Ne' 'Te' 'Th' 'song' 'mlnoise'};
CallCatLRI={ 'Ag' 'Be' 'DC' 'Di' 'LT' 'Ne' 'Te' 'Th' 'song'};
NCat11=length(CallCatRate);
NCat9=length(CallCatLRI);
NCat10=length(CallCatDP);

input_dir=pwd;
Subjects = dir(input_dir);

% Calculate the total number of units
Nunits=0;
for ss=1:length(Subjects)
    Indiv=Subjects(ss).name;
    if length(Indiv)==11
        Idir=fullfile(input_dir, Indiv);
        matfiles=dir(fullfile(Idir,'LRI_DP*.mat'));
        lm=length(matfiles);
        SS_Ind=zeros(lm,1);
        for ff=1:lm
            if ~isempty(strfind(matfiles(ff).name, 'ss'))
                SS_Ind(ff)=1;
            end
        end
        if strcmp(UT, 'SS')
            Indices=find(SS_Ind);
        elseif strcmp(UT, 'S')
            Indices=find(~SS_Ind);
        elseif strcmp(UT, 'ALL')
            Indices=1:lm;
        end
        LM=length(Indices);
        Nunits = LM + Nunits;
    end
end




% Initialize matrices for average dprimes, StdDprimes, Dsel, AvgRcallcat,
% StdRcallcat
MAvgDprime=NaN(Nunits, NCat10);%dprimes calculated between each category of vocalization and mlnoise with mlnoise
MStdDprime=NaN(Nunits, NCat10);
MAvgRcallcat=NaN(Nunits, NCat11);%calculated for as much sound categories as there are
MStdRcallcat=NaN(Nunits, NCat11);
UnitLRI=NaN(Nunits, NCat9);%calculated only for vocalization whine excepted
UnitSSI=NaN(Nunits, NCat9);
Unitnames=cell(Nunits,1);
Unitsubject=cell(Nunits,1);
UnitPath=cell(Nunits,1);


% Loop through units matfiles and fill in matrices
unit=0;
for ss=1:length(Subjects)
    Indiv= Subjects(ss).name;
    if length(Indiv)==11
        Idir=fullfile(input_dir, Indiv);
        calfiles=dir(fullfile(Idir,'LRI_DP*.mat'));
        lm=length(calfiles);
        SS_Ind=zeros(lm,1);
        for ff=1:lm
            if ~isempty(strfind(calfiles(ff).name, 'ss'))
                SS_Ind(ff)=1;
            end
        end
        if strcmp(UT, 'SS')
            Indices=find(SS_Ind);
        elseif strcmp(UT, 'S')
            Indices=find(~SS_Ind);
        elseif strcmp(UT, 'ALL')
            Indices=1:lm;
        end
        LM=length(Indices);
        fprintf('Loading Matfiles of %s\n',Indiv);
        for cc=1:LM
            unit = unit + 1;
            calfilename=calfiles(Indices(cc)).name;
            calfile=fullfile(Idir, calfilename);
            fprintf('Loading Matfile %s\n', calfiles(Indices(cc)).name);
            Cal=load(calfile);
            
            % checking the type of vocalization category we have here
            %first for Average spike rate and std spike rate
            outputRate=1:NCat11;
            inputRate=zeros(1,length(Cal.StimType));
            ii=0;
            for vv=1:NCat11
                ct = CallCatRate(vv);
                selector=find(strcmp(ct , Cal.StimType));
                if isempty(selector)
                    II=find(outputRate==vv);
                    outputRate(II)=[];
                elseif length(selector)==1
                    ii=ii+1;
                    inputRate(ii)=selector;
                else
                    fprintf('could not correctly attribute call type for average call rate\n');
                end
            end
            
            %second for SSI and LRI indices
            outputI=1:NCat9;
            inputI=zeros(1,length(Cal.StimTypeR));
            ii=0;
            for vv=1:NCat9
                ct = CallCatLRI(vv);
                selector=find(strcmp(ct , Cal.StimTypeR));
                if isempty(selector)
                    II=find(outputI==vv);
                    outputI(II)=[];
                elseif length(selector)==1
                    ii=ii+1;
                    inputI(ii)=selector;
                else
                    fprintf('could not correctly attribute call type for LRI and SSI\n');
                end
            end
            
            %finally for dprimes
            outputDP=1:NCat10;
            inputDP=zeros(1,length(unique(Cal.AllComptype)));
            ii=0;
            for vv=1:NCat10
                ct = CallCatDP(vv);
                selector=find(strcmp(ct , unique(Cal.AllComptype)));
                if isempty(selector)
                    II=find(outputDP==vv);
                    outputDP(II)=[];
                elseif length(selector)==1
                    ii=ii+1;
                    inputDP(ii)=selector;
                else
                    fprintf('could not correctly attribute call type for dprimes\n');
                end
            end
            
            
            % filling matrices with results
            MAvgRcallcat(unit,outputRate)= Cal.AvgRcallcat(inputRate);
            MStdRcallcat(unit,outputRate)= Cal.StdRcallcat(inputRate);
            UnitLRI(unit,outputI)=Cal.LRIcat(inputI);
            UnitSSI(unit,outputI)=Cal.SSIcat(inputI);
            MAvgDprime(unit, outputDP)=Cal.AvgDprime(inputDP);
            MStdDprime(unit, outputDP)=Cal.StdDprime(inputDP);
            Unitnames{unit}=calfilename;
            Unitsubject{unit}=Indiv;
            UnitPath{unit}=calfile;
            
        end
    
    fprintf('Done pulling data for %s\n',Subjects(ss).name);
    end
end

             
 % Save matrices in a matfile
 if ismac()
        [status username] = system('who am i');
        if strcmp(strtok(username), 'frederictheunissen')
            if strncmp('/auto/fdata/solveig',stim_name, 19)
            elseif strncmp('/auto/fdata/julie',stim_name, 17)
                Dfilename = fullfile('/Users','frederictheunissen','Documents','Data','Julie','matfile','AllUnits_LRI_DP.mat');
            end
        elseif strcmp(strtok(username), 'elie')
            Dfilename = fullfile('/Users','elie','Documents','MATLAB','data','matfile','AllUnits_LRI_DP.mat');
        end
else
    Dfilename=fullfile('/auto','k6','julie','calmatfile','AllUnits_LRI_DP.mat');
end
save(Dfilename, 'M*', 'U*', 'CallCatRate', 'CallCatDP', 'CallCatLRI');
fprintf('Done compiling results to %s\n', Dfilename);
return