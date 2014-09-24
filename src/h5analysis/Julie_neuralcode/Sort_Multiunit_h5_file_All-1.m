addpath(genpath('/auto/fdata/julie/Matlab'));
addpath(genpath('/auto/k2/share/tlab'));
addpath(genpath('/auto/k5/matlab714/toolbox/stats/stats'));
%cd /auto/fdata/julie/h5
cd /auto/k6/julie/h5
input_dir=pwd;
Subjects = dir(input_dir);
Subj=(3:8);
SignifSections=cell(6,1);
SignifSectionsDeg=cell(6,1);
NbSections=cell(6,1);
pTHRESHSections=cell(6,1);
for SS=1:6
    ss=Subj(SS);
    Idir=fullfile(input_dir, Subjects(ss).name);
    h5files=dir(fullfile(Idir,'*.h5'));
    %matfiles=dir(fullfile(Idir,'ZS*.mat'));
    LH=length(h5files);
    %LM=length(matfiles);
    SignifSections_local=zeros(LH,1);
    SignifSectionsDeg_local=zeros(LH,1);
    NbSections_local=zeros(LH,1);
    pTHRESHSections_local=zeros(LH,1);
    fprintf('Sorting h5 of %s\n',Subjects(ss).name);
    nh=0;
    for hh=1:LH
        h5file=fullfile(Idir, h5files(hh).name);
        fprintf('Calculating pValues of %s\n', h5files(hh).name);
        [nsections, pTHRESH, Section_pvalue, NSignif_Sections, pTHRESH_Deg, NSignif_Sections_Deg]=Sort_Multiunit_h5_file(h5file);
        if nsections ~= 0
            nh = nh+1;
            SignifSections_local(nh)=NSignif_Sections;
            SignifSectionsDeg_local(nh)=NSignif_Sections_Deg;
            NbSections_local(nh)=nsections;
            pTHRESHSections_local(nh)=pTHRESH;
        else
            fprintf('No value calculated for this file %s\n', h5files(hh).name);
        end
        fprintf('done\n');
    end
    SignifSections{SS}=SignifSections_local;
    SignifSectionsDeg{SS}=SignifSectionsDeg_local;
    NbSections{SS}=NbSections_local;
    pTHRESHSections{SS}=pTHRESHSections_local;
end
Res = struct();
Res.SignifSections=SignifSections;
Res.SignifSecionsDeg=SignifSectionsDeg;
Res.NbSections=NbSections;
Res.pTHRESHSections=pTHRESHSections;
if ismac()
    [status username] = system('who am i');
    if strcmp(strtok(username), 'frederictheunissen')
        if strncmp('/auto/fdata/solveig',stim_name, 19)
        elseif strncmp('/auto/fdata/julie',stim_name, 17)
            filename = fullfile('/Users','frederictheunissen','Documents','Data','Julie','matfile','Sort_Multiunith5_SignifSections.mat');
        end
    elseif strcmp(strtok(username), 'elie')
        filename = fullfile('/Users','elie','Documents','MATLAB','data','matfile','Sort_Multiunith5_SignifSections.mat');
    end
else
    filename=fullfile('/auto','k6','julie','matfile','Sort_Multiunith5_SignifSections2.mat');
end
save(filename, '-struct', 'Res');
fprintf('saved data under: %s\n', filename);
exit