%% sort the h5 according to the value of coherence information
addpath(genpath('/auto/fhome/julie/matlab/tlab'));
addpath(genpath('/auto/k5/matlab714/toolbox/stats/stats'));
%cd /auto/fdata/julie/h5
cd /auto/k6/julie/h5
input_dir=pwd;
Subjects = dir(input_dir);
Subj=(3:8);
fid_outF = fopen('h52Sort_FET.txt', 'wt');
fid_outM = fopen('h52Sort_MS.txt', 'wt');
fid_outT = fopen('h52Sort_TL.txt', 'wt');
fid_outH = fopen('h52Sort_HS.txt', 'wt');
fid_outJ = fopen('h52Sort_JEE.txt', 'wt');
for SS=1:6
    ss=Subj(SS);
    Idir=fullfile(input_dir, Subjects(ss).name);
    h5files=dir(fullfile(Idir,'*.h5'));
    LH=length(h5files);
    WhoList=repmat(randperm(5), 1, ceil(LH/5));
    mkdir(Idir, 'NonSignificant_h5');
    ih=0;
    for hh=1:LH
        nameh5=h5files(hh).name;
        h5file=fullfile(Idir, nameh5);
        fprintf('reading file %s\n', h5file);
        unit=read_unit_sectionsdata(h5file);
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
                if CoherenceInfo >= 3
                    ih=ih+1;
                    Who=WhoList(ih);
                    if Who==1
                        fprintf(fid_outF,'%s\n',h5file);
                    elseif Who==2
                        fprintf(fid_outM,'%s\n',h5file);
                    elseif Who==3
                        fprintf(fid_outT,'%s\n',h5file);
                    elseif Who==4
                        fprintf(fid_outH,'%s\n',h5file);
                    elseif Who==5
                        fprintf(fid_outJ,'%s\n',h5file);
                    end
                else
                    movefile (h5file, strcat(Idir, '/NonSignificant_h5/', nameh5));
                    fprintf(1,'the following h5 file is not significant and was moved in the NonSignificant_h5 folder: %s\n',h5file);
                end
            else
                movefile (h5file, strcat(Idir, '/NonSignificant_h5/', nameh5));
                fprintf(1,'the following h5 file has no coherence info data and was moved in the NonSignificant_h5 folder: %s\n',h5file);
            end
        else
            movefile (h5file, strcat(Idir, '/NonSignificant_h5/', nameh5));
            fprintf(1,'the following h5 file has no Call protocol data and was moved in the NonSignificant_h5 folder: %s\n',h5file);
        end
    end
end

fclose(fid_outF);
fclose(fid_outM);
fclose(fid_outT);
fclose(fid_outH);
fclose(fid_outJ);
exit
        