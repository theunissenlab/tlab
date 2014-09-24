function Sort_Multiunit_h5_file_all_SUBJ(SUBJ)
%% calculate the number of significant sections, the coherence and the...
...inter spike intervals for each h5 and store it...
...in the h5 and under Sort_Multiunith5_SignifSections.mat

addpath(genpath('/auto/fhome/julie/matlab/tlab'));
addpath(genpath('/auto/fhome/julie/matlab/STRFLab'));
addpath(genpath('/auto/k5/matlab714/toolbox/stats/stats'));
%cd /auto/fdata/julie/h5
cd /auto/k6/julie/h5
input_dir=pwd;

Subjects = dir(input_dir);

for ss=1:length(Subjects)
    Indiv=Subjects(ss).name;
    if strcmp(SUBJ, Indiv)
        Idir=fullfile(input_dir, Indiv);
        h5files=dir(fullfile(Idir,'Site*ss*.h5'));
        %matfiles=dir(fullfile(Idir,'ZS*.mat'));
        LH=length(h5files);
        %LM=length(matfiles);
        fprintf('Sorting h5 of %s\n',Subjects(ss).name);
        for hh=1:LH
            h5file=fullfile(Idir, h5files(hh).name);
            fprintf('Calculating pValues of %s\n', h5files(hh).name);
            [nsections, pTHRESH, Section_pvalue, NSignif_Sections, NSignif_Sections_rand, NSignif_Sections_FakeTrial, pTHRESH_Deg, NSignif_Sections_Deg, NSignif_Sections_Deg_rand, NSignif_Sections_Deg_FakeTrial]=Sort_Multiunit_h5_file(h5file);
            if nsections == 0
                fprintf('No value calculated for this file %s\n', h5files(hh).name);
            end
            fprintf('done\n');
            fprintf('Calculating coherence of %s\n', h5files(hh).name);
            compute_coherence_for_unit(h5file)
            fprintf('done\n');
            fprintf('Calculating spike time intervals of %s\n', h5files(hh).name);
            spike_auto_correlation_h5(h5file, 'Call', 'a');
            fprintf('done\n');
        end
    end
end

exit