function []=MetfilesScript_Indiv_JEE(SubjName)
addpath(genpath('/auto/fhome/julie/matlab/tlab'));
addpath(genpath('/auto/k5/matlab714/toolbox/stats/stats'));
%cd /auto/fdata/julie/h5
cd /auto/k6/julie/h5
input_dir=pwd;
Subjects = dir(input_dir);
for subj=1:length(Subjects)
    if strcmp(Subjects(subj).name, SubjName)
        ss=subj;
    end
end
%Subj=(3:8);
%for SS=1:6
    %ss=Subj(SS);
    Idir=fullfile(input_dir, Subjects(ss).name);
    h5files=dir(fullfile(Idir,'*.h5'));
    LH=length(h5files);
    %LM=length(matfiles);
    fprintf('Creating Matfile for %s\n',Subjects(ss).name);
    for hh=1:LH
        h5file=fullfile(Idir, h5files(hh).name);
        fprintf('Calculating Matfile of %s\n', h5files(hh).name);
        Matfile_construct_JEE(h5file);
        fprintf('done\n');
        %MatfilePath=fullfile(Idir, matfiles(hh).name);
        %fprintf('Calculating dprimes for %s\n', MatfilePath);
        %calls_selectivity(MatfilePath)
        %fprintf('done\n');
        clear h5file
    end
    
addpath(genpath('/auto/k5/matlab714/toolbox/stats/classreg'));
addpath(genpath('/auto/k5/matlab714/toolbox/econ'));
cd /auto/k6/julie/matfile
input_dir=pwd;
Idir=fullfile(input_dir, Subjects(ss).name);
matfiles=dir(fullfile(Idir,'ZS*.mat'));
LM=length(matfiles);
fprintf('Creating Calfile for %s\n',Subjects(ss).name);
for hh=1:LM
    MatfilePath=fullfile(Idir, matfiles(hh).name);
    fprintf('Calculating dprimes, Confusion matrices, gamma values for %s\n', MatfilePath);
    calls_selectivity_Modelonly(MatfilePath)
    fprintf('done making calculus on %s\n', MatfilePath);
    clear MatfilePath
end
    
%end
exit
end
