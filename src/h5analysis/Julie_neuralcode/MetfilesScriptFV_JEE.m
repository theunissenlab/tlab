function MetfilesScriptFV_JEE(UT,OW)
addpath(genpath('/auto/fhome/julie/matlab/tlab'));
addpath(genpath('/auto/k5/matlab714/toolbox/stats/stats'));
%cd /auto/fdata/julie/h5
cd /auto/k6/julie/h5

%UT is the unit Type can be 'S', 'SS' or 'ALL'
%OW is the overwriting switch can be 1 (matfiles of units that were already
%processed will be overwritten) or 0

if nargin==0
    UT='ALL';
    OW=1;
end
if nargin==1
    OW=1;
end

input_dir=pwd;
output_dir='/auto/k6/julie/matfile';
Subjects = dir(input_dir);
ExistingFiles=0;

for ss=1:length(Subjects)
    Indiv=Subjects(ss).name;
    if length(Indiv)==11
        fprintf('Looking for sound file names in the Vocalization Bank for %s\n', Indiv);
        [OldWav]=relate_sound_files(Indiv);
        Idir=fullfile(input_dir, Indiv);
        Outdir=fullfile(output_dir, Indiv);
        matfiles=dir(fullfile(Outdir,'FirstVoc*.mat'));
        LM=length(matfiles);
        
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
        
        fprintf('Creating Matfile for %s\n',Indiv);
        for hh=1:LH
            h5name=h5filesall(Indices(hh)).name;
            h5file=fullfile(Idir, h5name);
            
            %Make sure the file does not already exists
            FE=0;
            for mm = 1:LM
                if strcmp(h5name(1:(end-3)), matfiles(mm).name(10:(end-4)));
                    fprintf('%s has already a Matfile\n', h5file);
                    if OW==0
                        FE=1;%switch to file exist only if OverWriting is requested
                    end
                    ExistingFiles=ExistingFiles+1;
                end
            end
            if FE==0
                fprintf('Calculating Matfile of %s\n', h5file);
                Matfile_construct_FV_JEE(h5file,OldWav);
                fprintf('done\n');
            end
            
        end
        
    end
end
if OW==1
    sprintf('%d matfiles already existed and has been overwritten', ExistingFiles);
elseif OW==0
    sprintf('%d matfiles already existed and were kept unchanged', ExistingFiles);
end
exit
