function CalMatfilesScript_JEE(UT,AT)

addpath(genpath('/auto/fhome/julie/matlab/tlab'));
addpath(genpath('/auto/k5/matlab714/toolbox/stats/stats'));
addpath(genpath('/auto/k5/matlab714/toolbox/stats/classreg'));
addpath(genpath('/auto/k5/matlab714/toolbox/econ'));
%cd /auto/fdata/julie/h5
cd /auto/k6/julie/matfile

if nargin == 0
    UT='ALL';
end
if nargin == 1
    AT = 'R'; % 'R' for rate analysis running calls_selectivity_LRI_dprime set to 'CM' to run calculus of confusion matrices
end

input_dir=pwd;
Subjects = dir(input_dir);

for ss=1:length(Subjects)
    Indiv=Subjects(ss).name;
    if length(Indiv)==11
        Idir=fullfile(input_dir, Indiv);
        if strcmp(AT, 'R')
            matfiles=dir(fullfile(Idir,'WholeVoc*.mat'));
        elseif strcmp(AT, 'CM')
            matfiles=dir(fullfile(Idir,'FirstVoc*.mat'));
        end
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
        
        fprintf('Creating Calfile for %s\n',Indiv);
        for hh=1:LM
            MatfilePath=fullfile(Idir, matfiles(Indices(hh)).name);
            if strcmp(AT, 'R')
            	fprintf('Calculating dprimes, Index of selectivity for %s\n', MatfilePath);
                [Calfilename]=calls_selectivity_LRI_dprime(MatfilePath);
            elseif strcmp(AT, 'CM')
                fprintf('Calculating confusion matrices for %s\n', MatfilePath);
                [Calfilename]=calls_selectivity_ConfusionMatrix(MatfilePath);
            end
            
            fprintf('done making calculus on %s\nData save under %s\n', MatfilePath, Calfilename);
            clear MatfilePath
        end
    end
end
exit
