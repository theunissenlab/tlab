cd /auto/k6/julie/matfile
resultsDirectory='/auto/k6/julie/matfile';
addpath('/auto/k1/queued');
OW=0; %set Overwriting to 1 if you want to overwrite existing output files

input_dir=pwd;
Subjects = dir(input_dir);
cmd = 'calls_selectivity_ConfusionMatrix(''%s'');';
ExistingFiles=0;

for ss=1:length(Subjects)
    Indiv=Subjects(ss).name;
    if length(Indiv)==11
%     if strcmp(Indiv, 'GreBlu9508M')
        Idir=fullfile(input_dir, Indiv)
        matfiles=dir(fullfile(Idir,'FirstVoc*.mat'));
        outputmatfiles=dir(fullfile(Idir,'ConfMat_*.mat'));
        LO=length(outputmatfiles);
%        matfiles=dir(fullfile(Idir,'FirstVoc*Site2*ss1*.mat'));
        lm=length(matfiles);
        SS_Ind=zeros(lm,1);
        for ff=1:lm
            if ~isempty(strfind(matfiles(ff).name, 'ss'))
                SS_Ind(ff)=1;
            end
            %if ~isempty(strfind(matfiles(ff).name, 'e12'))
%                 SS_Ind(ff)=1;
%             end
%             if ~isempty(strfind(matfiles(ff).name, 'e21'))
%                 SS_Ind(ff)=1;
%             end
        end
        Indices=find(SS_Ind);
        LM=length(Indices);
        %Randlist=randperm(LM);
        for hh=1:LM
            MATName=matfiles(Indices(hh)).name;
            FE=0;
            for mm = 1:LO
                if strcmp(MATName(10:(end-4)), outputmatfiles(mm).name(9:(end-4)));
                    fprintf('%s has already a Matfile\n', MATName);
                    if OW==0
                        FE=1;%switch to file exist only if OverWriting is not requested
                    end
                    ExistingFiles=ExistingFiles+1;
                end
            end
            if FE==0
                jobParams = struct;
                jobParams.partition = 'all';
                jobParams.cpus = 2;
                jobParams.memory = 500;
                jobParams.out = fullfile(resultsDirectory, Indiv, sprintf('slurm_out_%s.txt', MATName));
                jobParams.err = jobParams.out;
                MatfilePath=fullfile(Idir, MATName);
                icmd = sprintf(cmd, MatfilePath); 
                fprintf('%s. Calling slurm_sbatch with command %s\n',Indiv, icmd);
                slurm_sbatch(icmd,jobParams);
            end
            
        end
    end
end
exit
