addpath(genpath('/auto/fhome/julie/matlab/tlab'));
addpath(genpath('/auto/k5/matlab714/toolbox/stats/stats'));
%cd /auto/fdata/julie/h5
cd /auto/k6/julie/h5

input_dir=pwd;

% Rescuing files from WhiWhi4522M that were wrongly placed here because of
% an error path
UT='S';
fid_out=fopen('List_h5files_tosort_WhiWhi4522M.txt', 'wt');

Indiv = 'WhiWhi4522M';
Idir=fullfile(input_dir, Indiv, 'NonSignificant_h5');
h5filesall=dir(fullfile(Idir,'Site1*s*.h5'));
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
                       % 
                        fprintf(1,'not significant\n');
                    else
                        CoherenceInfo
                        movefile (h5file, strcat(Idir, nameh5));
                         fprintf(1,'file moved back to main directory\n');
                         % create the list of files to be sorted
                         fprintf(fid_out, '%s\n', h5file);
                         
                    end
                else
                    %movefile (h5file, strcat(Idir, '/NonSignificant_h5/', nameh5));
                    fprintf(1,'no coherence info data\n');
                end
            else
                %movefile (h5file, strcat(Idir, '/NonSignificant_h5/', nameh5));
                fprintf(1,'no Call protocol data\n');
            end
        end
        fclose(fid_out);


%% Rescue ss files (65) that have value of info between 2.3 and 3
UT='SS';
Subjects = dir(input_dir);

for ss=1:length(Subjects)
    Indiv=Subjects(ss).name;
    if length(Indiv)==11
        Idir=fullfile(input_dir, Indiv, 'NonSignificant_h5');
        
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
                       % movefile (h5file, strcat(Idir, '/NonSignificant_h5/', nameh5));
                        fprintf(1,'not significant\n');
                    else
                        CoherenceInfo
                        movefile (h5file, fullfile(input_dir, Indiv, nameh5));
                         fprintf(1,'file moved back to main directory\n');
                         
                    end
                else
                    %movefile (h5file, strcat(Idir, '/NonSignificant_h5/', nameh5));
                    fprintf(1,'no coherence info data\n');
                end
            else
                %movefile (h5file, strcat(Idir, '/NonSignificant_h5/', nameh5));
                fprintf(1,'no Call protocol data\n');
            end
        end
    end
end

%% Rename these 65 files :-( because I used strcat instead of fullfile line 141

Subjects = dir(input_dir);

for ss=1:length(Subjects)
    Indiv=Subjects(ss).name;
    if length(Indiv)==11
        Idir=fullfile(input_dir, Indiv);
        h5filesall=dir(fullfile(Idir,'NonSignificant_h5Site*s*.h5'));
        lh=length(h5filesall);
        for hh=1:lh
            nameh5=h5filesall(hh).name;
            h5file=fullfile(Idir, nameh5);
            movefile (h5file, fullfile(Idir, nameh5(18:end)));
            fprintf(1,'file %s renamed %s\n', nameh5, nameh5(18:end));
        end
    end
end
                          
            
        