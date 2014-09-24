
%cd /Users/elie/Documents/MATLAB/data/calmatfile
cd /auto/k6/julie/calmatfile
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
    Indiv=Subjects(ss).name;
    if length(Indiv)==11
        Idir=fullfile(input_dir, Indiv);
        calfiles=dir(fullfile(Idir,'Cal*_nc.mat'));
        LC=length(calfiles);
        Nunits = LC + Nunits;
    end
end

% Initialize matrices for average dprimes, StdDprimes, Dsel, AvgRcallcat,
% StdRcallcat
MDsel=NaN(Nunits,Ncomp);
MAvgDprime=NaN(Nunits, Ncomp);
MStdDprime=NaN(Nunits, Ncomp);
MAvgRcallcat=NaN(Nunits, NCat);
MStdRcallcat=NaN(Nunits, NCat);
UnitLRI=NaN(Nunits, NCat);
UnitSSI=NaN(Nunits, NCat);
UnitDprimeMeans=NaN(Nunits,1);
UnitDprimeMean_wi=NaN(Nunits,1);
UnitDprimeStd_wi= NaN(Nunits,1);
UnitDprimeMean_be= NaN(Nunits,1);
UnitDprimeStd_be= NaN(Nunits,1);

Unitnames=cell(Nunits);
UnitR2A=NaN(Nunits, 3);
UnitPValLRatio=NaN(Nunits, 2);
UnitPvalue=NaN(Nunits, 3);


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
                    fprintf('could not correctly attribute call type\n');
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
                selector=find(cell2mat(Cal.Dsel.comptype(:,2))==ct);
                if isempty(selector)
                    II=find(outputI_comp1==vv);
                    outputI_comp1(II)=[];
                elseif length(selector)==1
                    ii=ii+1;
                    inputI_comp1(ii)=selector;
                else
                    fprintf('could not correctly attribute comparison type\n');
                end
            end
            
            % filling matrix with results
            MDsel(unit,outputI_comp1) = Cal.Dsel.values(inputI_comp1);
            
            
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
                    fprintf('could not correctly attribute comparison type\n');
                end
            end
                
            
            MAvgDprime(unit,outputI_comp2)= cell2mat(Cal.AvgDprime(inputI_comp2));
            MStdDprime(unit,outputI_comp2)= cell2mat(Cal.StdDprime(inputI_comp2));
            
            UnitDprimeMeans(unit)=nanmean(MAvgDprime(unit,outputI_comp2));
            UnitDprimeMean_wi(unit)=Cal.Mean_DP_wi;
            UnitDprimeStd_wi= Cal.Std_DP_wi;
            UnitDprimeMean_be= Cal.Mean_DP_be;
            UnitDprimeStd_be= Cal.Std_DP_be;
            UnitR2A(unit, :)=Cal.R2A;
            UnitPValLRatio(unit, :)=Cal.PValLRatio(2:3);
            UnitPvalue(unit, :)=Cal.Pvalue;
            Unitnames{unit}=[Indiv calfiles(cc).name];
        end
    
    fprintf('Done pulling data for %s\n',Subjects(ss).name);
    end
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
    Dfilename=fullfile('/auto','k6','julie','calmatfile','AllUnits.mat');
end
save(Dfilename, 'M*', 'U*', 'Compcat', 'CallCat');
fprintf('Done compiling results to %s\n', Dfilename);
return