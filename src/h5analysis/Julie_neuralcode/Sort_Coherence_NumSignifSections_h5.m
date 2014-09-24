function Sort_Coherence_NumSignifSections_h5(UT)
%% go through all the h5 and extract the number of significant sections...
... + give the Coherence info for all the h5 and store under...
    ...InfoCoherence_NumSignifSections_h5.mat
if nargin==0
    UT='ALL';
end

addpath(genpath('/auto/fhome/julie/matlab/tlab'));
addpath(genpath('/auto/k5/matlab714/toolbox/stats/stats'));
%cd /auto/fdata/julie/h5
cd /auto/k6/julie/h5
input_dir=pwd;
Subjects = dir(input_dir);
%Subj=(3:8);
Select_h5=cell(6,1); %will contain coherence info on first column, then
...number of significant sections under Bonferroni correction in second
    ...column and number of significant sections under gradual Bonferroni correction
fid_out1 = fopen(sprintf('h5_listFiles_CallProtocol_%s.txt',UT), 'wt');
SS=0;
th=0;
for ss=1:length(Subjects)
    Indiv=Subjects(ss).name;
    if length(Indiv)==11
        SS=SS+1;
        Idir=fullfile(input_dir, Indiv);
        h5filesall=dir(fullfile(Idir,'Site*.h5'));
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
        Select_h5_local=zeros(LH,7);
        ih=0;
        for hh=1:LH
            th=th+1;
            ih=ih+1;
            h5file=fullfile(Idir, h5filesall(Indices(hh)).name);
            unit=read_unit_h5file(h5file, 'r');
            SortType=unit.sortType;
            fprintf('reading file %s\n', h5file);
            fprintf(fid_out1, '%d\t%d\t%s\t%s\n', th, ih,SortType, h5file);
            

            nclasses = length(unit.classes);
            classId = 0;
            for ic=1:nclasses
                prot=unit.classes{ic};
                [path, Name, ext]=fileparts(unit.source_directory);
                if strcmp('WhiWhi4522M', Name) && strcmp('Site1', unit.site)
                    if strcmp('Call1c', prot(1:6))
                        classId = ic;
                        break;
                    end
                elseif strcmp('Call', prot(1:4))
                    classId = ic;
                    break;
                end
            end
            if isfield(unit, 'extra_info')
                Select_h5_local(ih,2)=unit.extra_info.(prot).sections_data.nb_signif_sections_Bonf;
                Select_h5_local(ih,3)=unit.extra_info.(prot).sections_data.nb_signif_sections_BonfDeg;
                p_values = unit.extra_info.(prot).sections_data.sections_pval;
                Select_h5_local(ih,4)=sum(p_values(:,3)<0.01);
                Select_h5_local(ih,5)=unit.extra_info.(prot).inter_spike_interval.Spike_proba_1ms;
                Select_h5_local(ih,6)=unit.extra_info.(prot).inter_spike_interval.Spike_proba_1ms_theo;
                Select_h5_local(ih,7)=unit.extra_info.(prot).inter_spike_interval.Stim_rate;
                
                

            end
            if isfield(unit, 'class_info')
                Select_h5_local(ih,1)=unit.class_info.(prot).info;
            end
            
        end
        Select_h5{SS}=Select_h5_local;
    end
end
if ismac()
    [status username] = system('who am i');
    if strcmp(strtok(username), 'frederictheunissen')
        if strncmp('/auto/fdata/solveig',stim_name, 19)
        elseif strncmp('/auto/fdata/julie',stim_name, 17)
            filename = fullfile('/Users','frederictheunissen','Documents','Data','Julie','matfile',sprintf('InfoCoherence_NumSignifSections_h5_%s.mat', UT));
        end
    elseif strcmp(strtok(username), 'elie')
        filename = fullfile('/Users','elie','Documents','MATLAB','data','matfile',sprintf('InfoCoherence_NumSignifSections_h5_%s.mat', UT));
    end
else
    filename=fullfile('/auto','k6','julie','matfile',sprintf('InfoCoherence_NumSignifSections_h5_%s.mat', UT));
end
save(filename, 'Select_h5');
fprintf('saved data under: %s\n', filename);
fclose(fid_out1);
exit
        