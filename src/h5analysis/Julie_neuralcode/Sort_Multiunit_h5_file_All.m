function Sort_Multiunit_h5_file_All(UT)
%% calculate the number of significant sections, the coherence and the...
...inter spike intervals for each h5 and store it...
...in the h5 and under Sort_Multiunith5_SignifSections.mat

if nargin==0
    UT='ALL';
end
addpath(genpath('/auto/fhome/julie/matlab/tlab'));
addpath(genpath('/auto/fhome/julie/matlab/STRFLab'));
addpath(genpath('/auto/k5/matlab714/toolbox/stats/stats'));
%cd /auto/fdata/julie/h5
cd /auto/k6/julie/h5
input_dir=pwd;

Subjects = dir(input_dir);

SignifSections=cell(6,1);
SignifSectionsDeg=cell(6,1);
NbSections=cell(6,1);
pTHRESHSections=cell(6,1);
SignifSections_Rand=cell(6,1);
SignifSectionsDeg_Rand=cell(6,1);
SignifSections_FakeTrial=cell(6,1);
SignifSectionsDeg_FakeTrial=cell(6,1);
SS=0;
for ss=1:length(Subjects)
    Indiv=Subjects(ss).name;
    if length(Indiv)==11
        SS=SS+1;
        Idir=fullfile(input_dir, Indiv);
        h5filesall=dir(fullfile(Idir,'Site*s*.h5'));
        lh=length(h5filesall);
        SS_Ind=zeros(lh,1);
        for ff=1:lh
            if ~isempty(strfind(h5filesall(ff).name, 'ss'))
                SS_Ind(ff)=1;
            end
        end
        if strcmp(UT, 'SS')
            Indices=find(SS_Ind);
        elseif strcmp(UT, 'S')
            Indices=find(~SS_Ind);
        elseif strcmp(UT, 'ALL')
            Indices=1:lh;
        end
        LH=length(Indices);
        SignifSections_local=zeros(LH,1);
        SignifSectionsDeg_local=zeros(LH,1);
        NbSections_local=zeros(LH,1);
        pTHRESHSections_local=zeros(LH,1);
        SignifSections_Rand_local=zeros(LH,1);
        SignifSectionsDeg_Rand_local=zeros(LH,1);
        SignifSections_FakeTrial_local=zeros(LH,1);
        SignifSectionsDeg_FakeTrial_local=zeros(LH,1);
        fprintf('Sorting h5 of %s\n',Subjects(ss).name);
        nh=0;
        for hh=1:LH
            h5file=fullfile(Idir, h5filesall(Indices(hh)).name);
            fprintf('Calculating pValues of %s\n', h5file);
            [nsections, pTHRESH, Section_pvalue, NSignif_Sections, NSignif_Sections_rand, NSignif_Sections_FakeTrial, pTHRESH_Deg, NSignif_Sections_Deg, NSignif_Sections_Deg_rand, NSignif_Sections_Deg_FakeTrial]=Sort_Multiunit_h5_file(h5file);
            if nsections ~= 0
                nh = nh+1;
                SignifSections_local(nh)=NSignif_Sections;
                SignifSectionsDeg_local(nh)=NSignif_Sections_Deg;
                SignifSections_Rand_local(nh)=NSignif_Sections_rand;
                SignifSectionsDeg_Rand_local(nh)=NSignif_Sections_Deg_rand;
                SignifSections_FakeTrial_local(nh)=NSignif_Sections_FakeTrial;
                SignifSectionsDeg_FakeTrial_local(nh)=NSignif_Sections_Deg_FakeTrial;
                NbSections_local(nh)=nsections;
                pTHRESHSections_local(nh)=pTHRESH;
            else
                fprintf('No value calculated for this file %s\n', h5file);
            end
            fprintf('done\n');
            fprintf('Calculating coherence of %s\n', h5file);
            compute_coherence_for_unit(h5file)
            fprintf('done\n');
            fprintf('Calculating spike time intervals of %s\n', h5file);
            spike_auto_correlation_h5(h5file, 'Call', 'a');
            fprintf('done\n');
        end
        SignifSections{SS}=SignifSections_local;
        SignifSectionsDeg{SS}=SignifSectionsDeg_local;
        NbSections{SS}=NbSections_local;
        pTHRESHSections{SS}=pTHRESHSections_local;
        SignifSections_Rand{SS}=SignifSections_Rand_local;
        SignifSectionsDeg_Rand{SS}=SignifSectionsDeg_Rand_local;
        SignifSections_FakeTrial{SS}=SignifSections_FakeTrial_local;
        SignifSectionsDeg_FakeTrial{SS}=SignifSectionsDeg_FakeTrial_local;
    end
end
Res = struct();
Res.SignifSections=SignifSections;
Res.SignifSecionsDeg=SignifSectionsDeg;
Res.NbSections=NbSections;
Res.pTHRESHSections=pTHRESHSections;
Res.SignifSections_Rand=SignifSections_Rand;
Res.SignifSectionsDeg_Rand=SignifSectionsDeg_Rand;
Res.SignifSections_FakeTrial=SignifSections_FakeTrial;
Res.SignifSectionsDeg_FakeTrial=SignifSectionsDeg_FakeTrial;
if ismac()
    [status username] = system('who am i');
    if strcmp(strtok(username), 'frederictheunissen')
        if strncmp('/auto/fdata/solveig',stim_name, 19)
        elseif strncmp('/auto/fdata/julie',stim_name, 17)
            filename = fullfile('/Users','frederictheunissen','Documents','Data','Julie','matfile','Sort_Multiunith5_Auditory_SignifSections.mat');
        end
    elseif strcmp(strtok(username), 'elie')
        filename = fullfile('/Users','elie','Documents','MATLAB','data','matfile','Sort_Multiunith5_Auditory_SignifSections.mat');
    end
else
    filename=fullfile('/auto','k6','julie','matfile','Sort_Multiunith5_Auditory_SignifSections.mat');
end
save(filename, '-struct', 'Res');
fprintf('saved data under: %s\n', filename);
exit