%% This code walk through all h5 files on k6/julie to find out any problem...
... of spike sorting or number of presentation of stim on recording day
%cd /Users/elie/Documents/MATLAB/data/h5
addpath(genpath('/auto/fhome/julie/matlab/tlab'));
cd /auto/k6/julie/h5
input_dir=pwd;
Subjects = dir(input_dir);
%Subj=(1:length(Subjects));
Fidoutssneed = fopen(fullfile(input_dir, 'List_h5_files_tosort.txt'), 'wt');
FidoutAuthorError = fopen(fullfile(input_dir, 'List_h5_files_AuthorError.txt'), 'wt');
FidoutTrialError = fopen(fullfile(input_dir, 'List_h5_files_TrialError.txt'), 'wt');

for ss=1:length(Subjects)
    %ss=Subj(SS);
    Indiv=Subjects(ss).name;
    if length(Indiv)==11
        Fidoutname=['List_h5files_' Indiv '.txt'];
        Idir=fullfile(input_dir, Indiv);
        Fidout = fopen(fullfile(Idir,Fidoutname), 'wt');
        allh5files=dir(fullfile(Idir,'*.h5'));
        LH=length(allh5files);
    
        for hh=1:LH
            Ih5=allh5files(hh).name;
            if isempty(strfind(Ih5, 'ss'))
                ref=strfind(Ih5, '_');
                ref=ref(3)+2;
                Ename=Ih5(1:ref);
                h5_ssFamily=cell(5,2);
                iss=0;
                for hh_int=1:LH
                    Ih5_int=allh5files(hh_int).name;
                    if ~isempty(strfind(Ih5_int, Ename)) && ~isempty(strfind(Ih5_int, 'ss'))
                        iss=iss+1;
                        Ih5_int_tot=sprintf('%s/%s', Indiv, Ih5_int);
                        [a,b]=system(sprintf('ls -l %s',Ih5_int_tot));
                        IN=find(b==' ');
                        Author=b(IN(2)+1 : IN(3)-1);
                        fprintf(Fidout, '%s\t%s\t%d\n', Ih5, Ih5_int, iss);
                        h5_ssFamily{iss,1}=Ih5_int;
                        h5_ssFamily{iss,2}=Author;
                        %fprintf('Calculating coherence of %s\n', Ih5_int);
                        %h5file=fullfile(Idir, Ih5_int);
                        %compute_coherence_for_unit(h5file)
                    end
                end
                if iss==0
                    fprintf(Fidoutssneed, '%s\t%s\n', Indiv, Ih5);
                else %verify that all files have the same authorship
                    ER=0;
                    for in=1:iss
                        if ~strcmp(h5_ssFamily(1,2), h5_ssFamily(in,2))
                            ER=1;
                        end
                    end
                    if ER==1
                        fprintf(FidoutAuthorError, '%s\t%s\n', Indiv, Ih5);
                    end
                end
                
                %% verifiy that sound were presented 10 times
                unit = read_unit_h5file(fullfile(Idir,Ih5), 'r');

                % Select the protocol (SelX, Callx, Maskx, STRFx...)
                %number of different protocols run for that unit, identify if the one
                %you are asking for is existing in this unit and selecting this one in
                %responses
                nclasses = length(unit.classes);

                classId = 0;
                for ic=1:nclasses
                    prot=unit.classes{ic};
                    stimType='Call';
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
                
                weirdTrials=0;
                nineTrial=0;
                elevenTrial=0;
                
                
                if classId ~= 0
                    
                    responses = unit.class_responses.(prot);
                    % This is the number of sound files played
                    nfiles = length(responses);
                    for isound = 1:nfiles
                        response=responses{isound};
                        ntrials=length(response.trials);
                        if ntrials~=10
                            weirdTrials = weirdTrials + 1;
                            if ntrials==9
                                nineTrial=nineTrial+1;
                            elseif ntrials==11
                                elevenTrial=elevenTrial+1;
                            end
                            
                        end
                    end
                end
                if weirdTrials~=0
                    fprintf(FidoutTrialError,'%s\t%s\tnb weird trials: %d\t nb eleven: %d\t nb nine: %d\n', Indiv, Ih5, weirdTrials, elevenTrial, nineTrial);
                end
                        
                
            end
        end
        fclose(Fidout);
    end
end
fclose(Fidoutssneed);
fclose(FidoutAuthorError);
fclose(FidoutTrialError);



