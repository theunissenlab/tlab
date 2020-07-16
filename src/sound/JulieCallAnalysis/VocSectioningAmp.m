function VocSectioningAmp(Folder, pl)
%% This function cut sound extract to isolate syllables in vocalizations based on the enveloppe of the syllable. 
... The syllables are centered at the peak and padded to 200 ms.  

%Set the parameters for sound filtering 
songfilt_args.f_low = 250;
songfilt_args.f_high = 12000;
songfilt_args.db_att = 0;



% duration=0.2; % duration in s of cuts
% endMax = 0.6;  % find the maximum in the first 600 ms  - this is for the
% neuro.

duration=0.35; % duration in s of cuts
endMax = 4.0;  % find the maximum in the first 4 s

% TDTsamprate = 24414;
TDTsamprate = 44100;       % This is not the TDTsamprate but using it so that I don't have to make too many mods in the code.

% Generate filter
nfilt = 512;
call_filter = fir1(nfilt,[songfilt_args.f_low*2.0/TDTsamprate, songfilt_args.f_high*2.0/TDTsamprate]);
call_filter_short = fir1(64,[songfilt_args.f_low*2.0/TDTsamprate, songfilt_args.f_high*2.0/TDTsamprate]);

if nargin<2
    pl=0; %set to 0 for no graph; 1 if you want to see final results per stim; 
end

%% Find the sound files on the cluster or on a local mac machine.
if nargin<1

    if ismac()
        [status username] = system('who am i');
        if strcmp(strtok(username), 'frederictheunissen')
            DataBasePath='/Users/frederictheunissen/Google Drive/Data/Julie/FullVocalizationBank';
        elseif strcmp(strtok(username), 'elie')
            DataBasePath='/Users/elie/Documents/ManipBerkeley/Recordings/group recordings';
        end
    end
else
    DataBasePath=Folder;
end
SoundFiles=dir(fullfile(DataBasePath, '*.wav'));

%% Cutting based on minima and maxima of the enveloppe below and above a threshold 10%

% This is the number of sound files
nfiles = size(SoundFiles,1);

% Initialize some data
ncutsTot = 0;
soundCutsTot = cell(1,nfiles);
spectroCutsTot = cell(1,nfiles);
indCutsTot = cell(1,nfiles);
soundNameCuts = cell(1,nfiles);
birdNameCuts= cell(1,nfiles);
vocTypeCuts= cell(1,nfiles);

% Make a file to save power
powerFile = fullfile(DataBasePath, 'ampcheck2020.txt');
fidPower = fopen(powerFile, 'w');

%Loop through files, detect sound periods
for isound = 1:nfiles
    stim_name=SoundFiles(isound).name;
     fprintf(1, 'Processing sound file %d/%d: %s\n', isound, nfiles, stim_name);
     
    % Read sound file. 
    [sound_temp, samprate] = audioread(fullfile(DataBasePath, stim_name));
    
    % We are changing to the TDT samprate...
    if (samprate ~= TDTsamprate )
        sound_in = resample(sound_temp, TDTsamprate, samprate);
        samprate = TDTsamprate;  
    else
        sound_in = sound_temp;
    end
    
    % Read the callid.
    [SoundPath, SoundName, ext]=fileparts(stim_name);
    voctype=SoundName(19:20);
    birdname = SoundName(1:10);
    
    % Make stereo sounds - mono.
    if size(sound_in,2)==2
        sound_in=(sound_in(:,1) + sound_in(:,2))/2;
    end
    
    % Read the Bird info file
    fid = fopen('/Users/frederictheunissen/Google Drive/Data/Julie/FullVocalizationBank/Birds_List_Acoustic.txt', 'r');
    birdInfo = textscan(fid, '%s %s %s %s %s %d');
    nInfo = length(birdInfo{1});
    fclose(fid);
    
    % Find matching bird
    birdInfoInd = 0;
    for iinfo=1:nInfo
        if (strcmp(birdInfo{1}(iinfo), birdname) )
            birdInfoInd = iinfo;
            break;
        end
    end
    
    if birdInfoInd == 0
        fprintf(1, 'Warning: Could not find entry for bird %s. Skipping.\n', birdname);
        continue;
    end
    % Readjust amplitude
    recLevel = birdInfo{6}(birdInfoInd);
    corrValTheo = power(10, (double(recLevel) - 100)/50);  % See recLevelAnal.m
    sound_in = sound_in./corrValTheo;
    
    % Filter the sound - we use to call songfilt_call here - which makes
    % sense for the Neuro data.
    if length(sound_in) < 3*nfilt
        sound_in = filtfilt(call_filter_short, 1, sound_in);
    else
        sound_in = filtfilt(call_filter, 1, sound_in);
    end
    fprintf(fidPower, '%s %s %s %d %f\n', SoundName, birdname, voctype, recLevel, std(sound_in));
    
    % Find the enveloppe
    amp_env = enveloppe_estimator(sound_in, samprate, 20, samprate);
    indEndMax = min(length(amp_env), fix(endMax*samprate));
    max_amp = max(amp_env(1:indEndMax));
    
    % Find the maxima above maxMinTh and minima below
    max_min_ind = [];
    j = 1;
    max_min_ind(j) = -1;
    maxMinTh = 0.1;
    nt = length(amp_env);
    for i=2:nt-1
        if ( (amp_env(i) > maxMinTh*max_amp) && (amp_env(i-1) < amp_env(i)) && (amp_env(i) > amp_env(i+1)) )
            j = j + 1;
            max_min_ind(j) = i;
        elseif ( (amp_env(i) <= maxMinTh*max_amp) && (amp_env(i) < amp_env(i-1)) && (amp_env(i) <= amp_env(i+1)) )
            j = j + 1;
            max_min_ind(j) = -i;
        end
    end
    j=j+1;
    max_min_ind(j) = -nt;

    
    % Clean up the max/min by keeping more extreme ones in succesive pairs
    i = 1;
    while i < length(max_min_ind)
        if (max_min_ind(i) > 0) && (max_min_ind(i+1) > 0)
            if amp_env(max_min_ind(i)) > amp_env(max_min_ind(i+1))
                max_min_ind(i+1) = [];
            else
                max_min_ind(i) = [];
            end
        elseif (max_min_ind(i) < 0) && (max_min_ind(i+1) < 0)
            if amp_env(-max_min_ind(i)) > amp_env(-max_min_ind(i+1))
                max_min_ind(i) = [];
            else
                max_min_ind(i+1) = [];
            end
        else
            i=i+1;
        end
    end
    
    if pl
        plot_sound_selection_amp(voctype, sound_in, samprate, amp_env, max_min_ind);
    end
    
    [ncuts soundCuts spectroCuts to fo indCuts] = cut_sound_selection_amp(sound_in, samprate, amp_env, max_min_ind, duration, pl);
    
    
    for ic=(ncutsTot+1):(ncutsTot+ncuts)
        soundNameCuts{ic} = SoundName;
        birdNameCuts{ic} = birdname;
        vocTypeCuts{ic} = voctype;
    end
    
    ncutsTot = ncutsTot + ncuts;
    soundCutsTot{isound} = soundCuts;
    spectroCutsTot{isound} = spectroCuts; 
    indCutsTot{isound} = indCuts;
    
     if pl
        pause();
    end
end

soundCutsAll = vertcat(soundCutsTot{:});
spectroCutsAll = vertcat(spectroCutsTot{:}); 
indCutsAll = vertcat(indCutsTot{:});

fclose(fidPower);

   
%% Save data to mat file
save('/Users/frederictheunissen/Google Drive/Data/Julie/FullVocalizationBank/vocCuts2020.mat', 'ncutsTot', 'to', 'fo', 'ncutsTot', 'soundCutsAll', 'spectroCutsAll', 'soundNameCuts', 'indCutsAll', 'birdNameCuts', 'vocTypeCuts', 'samprate', '-V7.3');

     