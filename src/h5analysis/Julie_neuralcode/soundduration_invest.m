
%This scripte retrieve the duration of the first cut in each stim and
%create and save a histogramm table
addpath(genpath('/auto/fhome/julie/matlab/tlab'));
addpath(genpath('/auto/k5/matlab714/toolbox/stats/stats'));
addpath(genpath('/auto/k5/matlab714/toolbox/stats/classreg'));
addpath(genpath('/auto/k5/matlab714/toolbox/econ'));
cd /auto/k6/julie/matfile
input_dir=pwd;
Subjects = dir(input_dir);
NVocSet = 0;

for ss=1:length(Subjects)
    Indiv=Subjects(ss).name;
    if length(Indiv)==11
        Idir=fullfile(input_dir, Indiv);
        fprintf(1,'%s\n', Indiv);
        matfiles=dir(fullfile(Idir,'WholeVoc_*.mat'));
        lm=length(matfiles);
        
        % select only one file per site
        Site_Ind=zeros(lm,1);
        site=1;
        for ff=1:lm
            if ~isempty(strfind(matfiles(ff).name, sprintf('Site%d',site)))
                fprintf(1,'file %s has been selected for site %d\n', matfiles(ff).name, site);
                site=site+1;
                Site_Ind(ff)=1;
            end
        end
        FinalIndices=find(Site_Ind);
        
        LS = length(FinalIndices);
        NVocSet=NVocSet + LS;
    end
end
HistDur=zeros(NVocSet*150,1);
Vocatype=cell(NVocSet*150,1);
vv=0;
for ss=1:length(Subjects)
    Indiv=Subjects(ss).name;
    if length(Indiv)==11
        Idir=fullfile(input_dir, Indiv);
        fprintf('%s\n', Indiv);        
        matfiles=dir(fullfile(Idir,'WholeVoc_*.mat'));
        lm=length(matfiles);
        
        % select only one file per site
        Site_Ind=zeros(lm,1);
        site=1;
        for ff=1:lm
            if ~isempty(strfind(matfiles(ff).name, sprintf('Site%d',site)))
                fprintf(1,'file %s has been selected for site %d\n', matfiles(ff).name, site);
                site=site+1;
                Site_Ind(ff)=1;
            end
        end
        FinalIndices=find(Site_Ind);
        
        LS = length(FinalIndices);
        
        for hh=1:LS
            MatfilePath=fullfile(Idir, matfiles(FinalIndices(hh)).name);
            fprintf(1,'Retrieving duration of first vocalization of each stim for %s\n', MatfilePath);
            Res=load(MatfilePath);
            Dur=Res.VocDuration;
            VocType=Res.VocType;
            DD=length(Dur);
            for dd=1:DD
                vv=vv+1;
                HistDur(vv)=Dur{dd}(1);
                Vocatype{vv}=VocType{dd};
            end
            fprintf(1,'done making calculus on %s\n', MatfilePath);
            clear MatfilePath
        end
    end
end
HistDur=HistDur(1:vv);
VocaType=Vocatype(1:vv);
filename=fullfile('/auto','k6','julie','matfile','HistoDuration1Voc.mat');
save(filename, 'HistDur','VocaType');
exit
