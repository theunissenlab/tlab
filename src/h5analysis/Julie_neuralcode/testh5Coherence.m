addpath(genpath('/auto/fdata/julie/Matlab'));
addpath(genpath('/auto/fhome/julie/matlab/tlab'));
addpath(genpath('/auto/fhome/julie/matlab/STRFLab'));
addpath(genpath('/auto/k5/matlab714/toolbox/stats/stats'));
%cd /auto/fdata/julie/h5
cd /auto/k6/julie/h5
input_dir=pwd;
Subjects = dir(input_dir);
Subj=(3:8);

for SS=1:6
    ss=Subj(SS);
    Idir=fullfile(input_dir, Subjects(ss).name);
    h5files=dir(fullfile(Idir,'*.h5'));
    %matfiles=dir(fullfile(Idir,'ZS*.mat'));
    LH=length(h5files);
    for hh=1:LH
        h5file=fullfile(Idir, h5files(hh).name);
        fprintf('Calculating coherence of %s\n', h5files(hh).name);
        compute_coherence_for_unit(h5file)
    end
end

