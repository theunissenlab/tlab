function Sort_h5files_CoherenceInfoBased_SS(UT)
%% sort the h5 according to the value of coherence information
addpath(genpath('/auto/fhome/julie/matlab/tlab'));
addpath(genpath('/auto/k5/matlab714/toolbox/stats/stats'));
%cd /auto/fdata/julie/h5
cd /auto/k6/julie/h5

if nargin==0
    UT='ALL';
end

fid=fopen('SortingTroubleAuditiveUnits.txt','w');

input_dir=pwd;
Subjects = dir(input_dir);

for ss=1:length(Subjects)
    Indiv=Subjects(ss).name;
    if length(Indiv)==11
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
            InfoThresh=2.3;
        elseif strcmp(UT, 'S')
            Indices=find(~SS_Ind);
            InfoThresh=3;
        elseif strcmp(UT, 'ALL')
            Indices=1:lh;
            InfoThresh=input('What threshold value do you want for sorting all the units?\n3 is used for s-files (non-sorted)\n 2.3 is used for ss-files (sorted files)\n');
        end
        
        LH=length(Indices);
        if exist('NonSignificant_h5', 'dir')==7
        else
            mkdir(Idir, 'NonSignificant_h5');
        end
        for hh=1:LH
            nameh5=h5filesall(Indices(hh)).name;
            h5file=fullfile(Idir, nameh5);
            fprintf('reading file %s\n', h5file);
            unit=read_unit_h5file(h5file, 'r');
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
            if classId~=0
                if isfield(unit, 'class_info')
                    CoherenceInfo=unit.class_info.(prot).info;
                    
                    if CoherenceInfo < InfoThresh
                        movefile (h5file, strcat(Idir, '/NonSignificant_h5/', nameh5));
                        fprintf(1,'the following h5 file is not significant and was moved in the NonSignificant_h5 folder: %s\n',h5file);
                    end
                else
                    movefile (h5file, strcat(Idir, '/NonSignificant_h5/', nameh5));
                    fprintf(fid,'%s no coherence info data and was moved in the NonSignificant_h5 folder\n',h5file);
                end
            else
                movefile (h5file, strcat(Idir, '/NonSignificant_h5/', nameh5));
                fprintf(fid,'%s no Call protocol data and was moved in the NonSignificant_h5 folder\n',h5file);
            end
        end
    end
end

fclose(fid);
exit
        