function CalMatfilesScript_JEE_Ind(Indiv)
addpath(genpath('/auto/fhome/julie/matlab/tlab'));
addpath(genpath('/auto/k5/matlab714/toolbox/stats/stats'));
addpath(genpath('/auto/k5/matlab714/toolbox/stats/classreg'));
addpath(genpath('/auto/k5/matlab714/toolbox/econ'));
%cd /auto/fdata/julie/h5
cd /auto/k6/julie/matfile
input_dir=pwd;

    if length(Indiv)==11
        Idir=fullfile(input_dir, Indiv);
        matfiles=dir(fullfile(Idir,'ZS*.mat'));
        LM=length(matfiles);
        fprintf('Creating Calfile for %s\n',Indiv);
        for hh=1:LM
            MatfilePath=fullfile(Idir, matfiles(hh).name);
            fprintf('Calculating dprimes, Index of selectivity, Confusion matrices for %s\n', MatfilePath);
            calls_selectivity_noconfusion(MatfilePath)
            fprintf('done making calculus on %s\n', MatfilePath);
            clear MatfilePath
        end
    end
exit
