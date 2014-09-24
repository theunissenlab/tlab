function SpectroCorrMat
% This function creates for each site a matrix of correlation values of all
% the spectro played back for all birds and sites

cd /auto/k6/julie/matfile
input_dir = pwd;
Subjects = dir(input_dir);
for ss=1:length(Subjects)
    Indiv=Subjects(ss).name;
    if length(Indiv)==11
        allFiles = dir(fullfile(input_dir, Indiv,'FirstVoc*.mat'));
        NF = length(allFiles);
        for nsite=1:8
            % find a matfile for each site
            Matsite='';
            for nf = 1:NF
                File_n=allFiles(nf).name;
                if str2num(File_n(14))==nsite
                    Matsite = allFiles(nf).name;
                    break
                end
            end
            if isempty(Matsite)
                fprintf(1, 'No file could be find for %s Site %d\n', Indiv, nsite);
            else
                % Retrieve the Matfile that contains the spectrograms
                MAT = load(fullfile(input_dir, Indiv,Matsite));

                % reorganize the vector of stims so that all the stims of the same category
                % are together
                VocTypeList=MAT.VocType;
                TDTWavfiles = MAT.TDT_wavfiles;
                VocTypeSel=cell(size(VocTypeList,1),1);
                StimIndices= nan(size(VocTypeList,1), 1);
                TDTwav = cell(size(VocTypeList,1),1);
                StimType=unique(VocTypeList);
                RestrictAna=find(~strcmp(StimType, 'Wh'));% localize the stims we don't want to study here and get rid of them
                if sum(strcmp(StimType, 'mlnoise'))%make sure that there are some mlnoise for that file. 
                    RestrictAna(find(RestrictAna==find(strcmp(StimType, 'mlnoise'))))=[];
                end
                StimTypeR = StimType(RestrictAna);
                NstimTypeR=length(StimTypeR);
                cc = 0;
                for vt=1:NstimTypeR
                    st=StimTypeR(vt);
                    selector=strcmp(VocTypeList, st);
                    selectorInd=find(selector);
                    NselectorInd=length(selectorInd);
                    for cc_temp=1:NselectorInd
                        cc=cc+1;
                        StimIndices(cc)=selectorInd(cc_temp);
                        VocTypeSel{cc}=VocTypeList{selectorInd(cc_temp)};
                        TDTwav{cc}=TDTWavfiles{selectorInd(cc_temp)};
                    end
                end
                StimIndices = StimIndices(1:cc);
                VocTypeSel = VocTypeSel(1:cc);
                TDTwav = TDTwav(1:cc);
                
                % set up the output matrix
                NSTIMS=length(StimIndices);
                CORR = nan(NSTIMS, NSTIMS);

                % loop through stims
                for s1=1:NSTIMS
                    Spectro1 = MAT.Spectro{StimIndices(s1)};
                    for s2=s1:NSTIMS
                        %fprintf(1,'s1=%d s2=%d\n', s1, s2);
                        Spectro2 = MAT.Spectro{StimIndices(s2)};
                        if ~length(Spectro1)==length(Spectro2)
                            fprintf(1,'WARNING, stims of different size %d and %d points\n', length(Spectro1), length(Spectro2));
                        end
                        %calculate amplitude= f(t) and find the best allignment ...
                        %between stims using xcorr(X,Y,'unbiased')
                        Spec1=reshape(Spectro1, length(MAT.Spectrofo), length(MAT.Spectroto));
                        Spec2=reshape(Spectro2, length(MAT.Spectrofo), length(MAT.Spectroto));
                        Amp1=sum(Spec1);
                        Amp2=sum(Spec2);
                        [C,LAGS] = xcorr(Amp1, Amp2, 'unbiased');
                        LAG = LAGS(find(C==max(C)));
                        LAG = LAG(1);
                        if LAG>=0
                            RfSpec1 = Spec1(:, LAG+1:end);
                            RfSpec2 = Spec2(:, 1:end-LAG);
                        else
                            RfSpec1 = Spec1(:, 1:end+LAG);
                            RfSpec2 = Spec2(:, 1-LAG:end);
                        end
                        CORR(s1, s2)= corr2(RfSpec1, RfSpec2);
                        if s1~=s2
                            CORR(s2, s1)=CORR(s1, s2);
                        end
                    end
                end

                % save the output matrix
                Res.CORR=CORR;
                Res.VocTypeSel=VocTypeSel;
                Res.TDTwav=TDTwav;
                Res.StimIndices= StimIndices;
                if ismac()
                    [status username] = system('who am i');
                    if strcmp(strtok(username), 'frederictheunissen')
                        if strncmp('/auto/fdata/solveig',stim_name, 19)
                        elseif strncmp('/auto/fdata/julie',stim_name, 17)
                            filename = fullfile('/Users','frederictheunissen','Documents','Data','Julie','matfile',Res.subject,['FirstVoc_' Res.Site '.mat']);
                        end
                    elseif strcmp(strtok(username), 'elie')
                        filename = fullfile('/Users','elie','Documents','MATLAB','data','matfile',Res.subject,['FirstVoc_' Res.Site '.mat']);
                    end
                else
                    filename=fullfile('/auto','k6','julie','matfile',Indiv,sprintf('CorrFirstVoc_Site%d.mat',nsite));
                end
                save(filename, '-struct', 'Res');
                fprintf('saved data under: %s\n', filename);
            end
        end
    end
end
                    
                    