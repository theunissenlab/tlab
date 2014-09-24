function [OldWav]=relate_sound_files(Indiv)

% Find on which machine we are and where are stored the data
    DataDir = '/auto/fdata/julie';
    if ismac()
            [status username] = system('who am i');
            if strcmp(strtok(username), 'frederictheunissen')
                if strcmp('/auto/fdata/julie',DataDir)
                    DataDir='/Users/frederictheunissen/Documents/Data/Julie';
                end
            elseif strcmp(strtok(username), 'elie')
                if strcmp('/auto/fdata/julie',DataDir)
                    DataDir='/Users/elie/Documents/MATLAB/data';
                    
                end
            end
    end
    
    % Find the protocol file that list the original wav files used as stims
    ProtoDir = strcat(DataDir, '/Stims/', Indiv);
    %ProtoDir = strcat('/auto/k8/fdata/julie/stim', Indiv);
    TrackDir=strcat(ProtoDir, '/Call');
    TrackFile=dir(fullfile(TrackDir, 'WavTracker*.txt'));
    FidTrack=fopen(fullfile(TrackDir, TrackFile.name));
    if FidTrack==-1
        fprintf('No readable Wavtracker for subject %s\n', Indiv);
        OldWav={};
    else
        header = textscan(FidTrack, '%s', 33);
        Trackdata = textscan(FidTrack, '%s', 'delimiter','\n', 'MultipleDelimsAsOne', 1);
        fclose(FidTrack);
        Nstims=size(Trackdata{1},1);
        OldWav=cell(Nstims, 5);
        ProtFile=dir(fullfile(ProtoDir, 'stimProtocol.txt'));
        FidProt = fopen(fullfile(ProtoDir, ProtFile.name));
        ProtData=textscan(FidProt, '%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s', 'delimiter', '\t');
        fclose(FidProt);
        nTDTstims=size(ProtData{1},1);
        for ii=1:Nstims
            Lineii=textscan(Trackdata{1}{ii}, '%s');
            len=size(Lineii{1},1);
            TDTname=Lineii{1}{1};
                stim_Idx=0;
            for jj=1:nTDTstims
                Protoname=ProtData{2}{jj};
                ProtonameParts = textscan(Protoname, '%s', 'Delimiter', '\\/');
                if strcmp(TDTname, ProtonameParts{1}{end})
                    stim_Idx=jj;
                end
            end

            if len==3
                OldWav{ii,1}=stim_Idx;
                OldWav{ii,2}=Lineii{1}{1};
                OldWav{ii,3}=Lineii{1}{3};
            elseif len==8
                OldWav{ii,1}=stim_Idx;
                OldWav{ii,2}=Lineii{1}{1};
                OldWav{ii,3}=Lineii{1}{3};
                OldWav{ii,4}=Lineii{1}{5};
                OldWav{ii,5}=Lineii{1}{7};
            else
                fprintf('Problem in the wavtracker of that subject with the stim %s, a space sign must be hiding somewhere\n', Lineii{1}{1});
            end
        end
    end
    
    
    
    