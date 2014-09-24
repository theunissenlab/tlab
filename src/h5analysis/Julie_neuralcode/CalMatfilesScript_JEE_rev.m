addpath(genpath('/auto/fhome/julie/matlab/tlab'));
addpath(genpath('/auto/k5/matlab714/toolbox/stats/stats'));
addpath(genpath('/auto/k5/matlab714/toolbox/stats/classreg'));
addpath(genpath('/auto/k5/matlab714/toolbox/econ'));
%cd /auto/fdata/julie/h5
cd /auto/k6/julie/matfile
input_dir=pwd;
Subjects = dir(input_dir);
Subj=[10 9 8 7 6 5 4 3 2 1];
for SS=1:length(Subjects)
    ss=Subj(SS);
    Indiv=Subjects(ss).name;
    if length(Indiv)==11
        Idir=fullfile(input_dir, Subjects(ss).name);
        matfiles=dir(fullfile(Idir,'ZS*.mat'));
        LM=length(matfiles);
        fprintf('Creating Calfile for %s\n',Subjects(ss).name);
        for hh=1:LM
            MatfilePath=fullfile(Idir, matfiles(hh).name);
            fprintf('Calculating dprimes, Index of selectivity, Confusion matrices for %s\n', MatfilePath);
            calls_selectivity_noconfusion(MatfilePath)
            fprintf('done making calculus on %s\n', MatfilePath);
            clear MatfilePath
        end
    end
end
AllUnits_calculator
exit
