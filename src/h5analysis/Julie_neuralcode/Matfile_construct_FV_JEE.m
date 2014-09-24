function [filename]=Matfile_construct_FV_JEE(h5Path, OldWav, pl)
%% This function cut tdt stim to isolate vocalizations as much as possible...
...and just take the first vocalization of each stim that it places in ...
    ...a window of xxxms

duration=0.02; % duration in ms of min silence period
if nargin<3
    pl=0; %set to 0 for no graph; 1 if you want to see final results per stim; 2 if you to see both steps graph
end


%%  I choose 600ms as a window to put the first voc of each stim in according...
... to values obtain with soundduration_invest.m
Win=0.6;

    %% Read input data from data base
unit = read_unit_h5file(h5Path, 'r');

%% Select the protocol (SelX, Callx, Maskx, STRFx...)
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

if classId ~= 0

    responses = unit.class_responses.(prot);

    %% Cutting procedure based on zero values to calculate intensity values of each cut
    % This is the number of sound files played
    nfiles = length(responses);
    
    % set up output cell arrays
    intensity = cell(nfiles,1);
    dur = cell(nfiles,1);
    voctype = cell(nfiles,1);
    Sound = cell(nfiles,1);
    histDur=[];
    
    % Find on which machine we are and where are stored the data
    [SoundPath] =fileparts(responses{1}.original_wavfile);
    SoundPathParts = textscan(SoundPath, '%s', 'Delimiter', '\\/');
    DataDir = '';
    k = 2;
    while ~strcmp(SoundPathParts{1}{k},'Stims')
        DataDir = strcat(DataDir, '/', SoundPathParts{1}{k});
        k=k+1;
    end
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
    
    %Loop through files, detect sound periods
    for isound = 1:nfiles
        response=responses{isound};
        stim_name=response.tdt_wavfile;
        stim_name = strcat(DataDir, stim_name(18:end));
        [sound_in, samprate] = wavread(stim_name);
        iintensity=[];
        idur=[];
        iSound = zeros(0);
    
       % Read the callid.
        if strcmp(response.stim_type, 'song')
            vocid=response.stim_type;
        elseif strcmp(response.stim_type, 'call')
            vocid=response.callid;
        else
            vocid=response.stim_type;
        end
    
        %if strcmp(callId,'Te')||strcmp(callId,'Th')||strcmp(callId,'DC')||strcmp(callId,'LT')
            % Divide the sound into areas of non-zero time
            ind = find(sound_in);
            on = zeros(size(sound_in));
            sound_p1 = zeros(size(sound_in));
            sound_p1(1:end-1) = sound_in(2:end);
            sound_m1 = zeros(size(sound_in));
            sound_m1(2:end) = sound_in(1:end-1);
            sound_diff = sound_p1-sound_m1;
            ind2 = find(sound_diff);
            
            % Sounds is "on" when value is not zero and derivative is not zero
            on(ind) = 0.08;    % It is set at 0.8 to make a nice graphic.
            on(ind2) = 0.08;
            if pl==2        
                plot_sound_spike_selection_h5(response, sound_in, samprate, on);
                pause(1);
            end
    
            % Find beggining and end
            laston = 0;
            for i=1:length(on)
                if (laston == 0) && (on(i)~=0)
                    begInd = i;
                elseif (laston~=0) && (on(i) == 0) || (laston~=0) && (i==length(on))
                    endInd = i;
                    durInd = endInd - begInd;
                
                    
                    if (durInd > samprate*duration)    % Only examine stimulus intervals of certain length
                    % calculate duration, peak frequency and intensity
                    durVal = (endInd-begInd)/samprate;
                    intVal = std(sound_in(begInd:endInd));
                    %[Pxx,f] = pwelch(sound_in(begInd:endInd), nw, fix(nw/2) , nw, samprate);
                    %Pmax = max(Pxx);
                    %fVal = f(Pxx==Pmax);
            
                     % Stuff values
                    iintensity = [iintensity intVal];
                    histDur = [histDur durVal];
                    idur = [idur durVal];
                    iSound = [iSound [begInd ; endInd]];
                    % fprintf(1, 'Beg= %.2f End= %.2f Int= %.2f Dur= %.2f Freq= %.0f RD= %.2f\n', begTime, endTime, 20*log10(intVal), durVal, fVal, rateDiffVal);
                    end
                end
                laston = on(i);
            end
            voctype{isound} = vocid;       
            intensity{isound}=iintensity;
            dur{isound}=idur;
            Sound{isound}=iSound;
    end
    figure(2)
    hist(histDur,60)



    %%  Define the length of the section
    Lsection=ceil(Win*samprate);



    %% Now run a second cutting system based on the intensity of the signal to get nicer cuts and enlarged the saved window around the vocalization
    % Create the structure that will contain the results
    Res = struct();
    %'subject',{},'Site', {},'VocType', {}, 'tdt_wavfiles', {},'original_wavfiles',{}, 'EndBegIndices', {}, 'Cut_orders', {}, 'sex', {}, 'age', {}, 'related', {}, 'trials', {}, 'PSTH',{}, 'spectro', {}, 'section_cat', {});

    % indicate info on the site
    [pathh5, nameh5, ext]=fileparts(h5Path);
    Res.subject=pathh5(end-10:end);
    Res.Site=nameh5;

    % prepare cell arrays
    coeff=8; %estimate of 8 sections per stim on average
    VocType=cell(nfiles*coeff,1);
    TDT_wavfiles=cell(nfiles*coeff,1);
    VocBank_Wavfiles=cell(nfiles*coeff,1);
    Original_wavfiles=cell(nfiles*coeff,1);
    SectionWave=cell(nfiles*coeff,1);
    WavIndices=cell(nfiles*coeff,1);
    SectionLength = zeros(nfiles*coeff,1);
    ESex=cell(nfiles*coeff,1);
    Eage=cell(nfiles*coeff,1);
    Erelated=cell(nfiles*coeff,1);
    Trials=cell(nfiles*coeff,1);
    Trials_BG=cell(nfiles*coeff,1);
    PSTH=cell(nfiles*coeff,1);
    PSTH_BG=cell(nfiles*coeff,1);
    Spectro=cell(nfiles*coeff,1);
    Section_cat=cell(nfiles*coeff,1);
    histDur2=[];
    n1section = 0;
    

    for isound = 1:nfiles
        response=responses{isound};
        stim_name=response.tdt_wavfile;
        stim_number=str2double(response.number);
        if ~isempty(OldWav)
            VocBank_idx= find(cell2mat(OldWav(:,1))==stim_number);
        end
    
        % Read the stim wave files on the cluster on a local mac machine.
        if ismac()
            [status username] = system('who am i');
            if strcmp(strtok(username), 'frederictheunissen')
                if strncmp('/auto/fdata/solveig',stim_name, 19)
                    stim_name = strcat('/Users/frederictheunissen/Documents/Data/solveig', stim_name(20:end));
                elseif strncmp('/auto/fdata/julie',stim_name, 17)
                    stim_name = strcat('/Users/frederictheunissen/Documents/Data/Julie', stim_name(18:end));
                end
            elseif strcmp(strtok(username), 'elie')
                if strncmp('/auto/fdata/solveig',stim_name, 19)
                    stim_name = strcat('/Users/frederictheunissen/Documents/Data/solveig', stim_name(20:end));
                elseif strncmp('/auto/fdata/julie',stim_name, 17)
                    stim_name = strcat('/Users/elie/Documents/MATLAB/data', stim_name(18:end));
                end
            end
        end
        [sound_in, samprate] = wavread(stim_name);
    
        %Find Silences longer than 60ms between vocalizations
        idur2=[];
        iintensity=intensity{isound};
        %%%%%%%%%%%%%%%Increase this a little more??? to enable separation
        %%%%%%%%%%%%%%%of nest calls?
        Thresh=0.35*mean(iintensity); %silence is defined as 35% of average intensity calculated above
        yy = 1;
        Duration = ceil(0.01*samprate);% calls have to be spaced by at least 10ms to be separated. 10ms tail is added after the end sensus stricto of each sound and before the begining
        Sil = zeros(0); %Initialize the matrix that will contain start and end indices of silence slots
        while (yy>0) && (yy<length(sound_in))
            if abs(sound_in(yy))<Thresh
                Beg = yy;
                while abs(sound_in(yy))<Thresh && (yy<length(sound_in))
                    yy = yy + 1;
                end
                End = yy;
                if End-Beg>(Duration)
                    Sil = [Sil [Beg ; End]];
                end
            end
            yy = yy + 1;
        end
        [r,c]=size(Sil);
        on=ones(size(sound_in));
        if c~=0
            for cc=1:c
                on((Sil(1,cc)+Duration):1:(Sil(2,cc)-Duration))=0;% indicate zones of silences (0) and zones around calls that should be kept (60ms after the end and 10ms before the begining)
            end
            if Sil(2,c)==length(on)%if the wavfile ends by a silence
                on((Sil(2,c)-Duration):1:end)=0;% then make sure the silence was not supress at the end by the previous statement
            end
        end
        if on(end)==1
            on= [on;ones(Duration,1)];%in case the stim end by a vocalization, add some space after (10ms)
        end
        on=on*0.8;
        if pl>0        
            plot_sound_spike_selection_h5(response, sound_in, samprate, on);
            pause;
        end
    
        %% Isolate spikes that relate to before the 1rst vocalization and calculate avg background spike rate
        
        ntrials=length(response.trials);
        bgdur=zeros(1,ntrials);
        spike_count_bg = zeros(1,ntrials);
        ntimebins = round(1000*Win); % Number of time bins in ms of the background psth window
        psthbg = zeros(1, ntimebins);
        BG_Trials_temp = cell(ntrials,1);
        for it=1:ntrials
            trial = response.trials{it};
            bgdur(it)=trial.pre_time;
            spike_time_trials = trial.spikeTimes;
            
            %calculate spike rate for he whoe background period
            spike_time1=(spike_time_trials<0);
            spike_time2=(spike_time_trials>=(-bgdur(it)));
            SpikesIndices=find(spike_time1 .* spike_time2);
            spike_count_bg(it)=length(SpikesIndices);
            %bg_mean_rate = mean(spike_count_bg)./(bgdur);
            %bg_std_rate = std(spike_count_bg)./(bgdur);
            spike_rate_bg = spike_count_bg./bgdur;
            
            %calculate a psth background for a 600ms window (Win) before
            %the stim
            spike_time3=(spike_time_trials>=(-Win));
            SpikesIndicesPSTHbg=find(spike_time1 .* spike_time3);
            Spike_time_PSTHbg=spike_time_trials(SpikesIndicesPSTHbg);
            BG_Trials_temp {it} = (Spike_time_PSTHbg - (-Win))*1000; % here spike times arrival are referred to the begining of the 600ms window instead of the begining of the stim and in ms instead of s
            ns=length(Spike_time_PSTHbg);
            spike_array = zeros(1, ntimebins);
            for is=1:ns
                time_ind = ceil(Spike_time_PSTHbg(is)*1000-(-Win*1000));
                if (time_ind < 1 || time_ind > ntimebins)
                    fprintf(1, 'Warning time index out of bounds for stim# %s trial %d: time_ind = %d ntimebins = %d\n', response.number, it, time_ind, ntimebins);
                    continue;
                end
                spike_array(time_ind) = spike_array(time_ind) +1;
            end
            psthbg = psthbg + spike_array;
        end
        psthbg = psthbg./ntrials;
        
        
    
        %% Find beggining and end of two first vocalization-section
        WavIndices_temp=cell(2,1);
        nsections = 0;
        laston = 0;
        begInd=0;
        i=0;
        while nsections<2 && i<length(on)
            nsections<2 && i<length(on);
            i = i + 1;
            %for i=1:length(on)
                if (laston == 0) && (on(i)~=0)
                    begInd = i;
                    
                elseif (laston~=0) && (on(i) == 0) || (laston~=0) && (i==length(on))
                    endInd = i;
                    durInd = endInd - begInd;
                    if begInd==1 %in case the stim start by a vocalization, the 10ms before begining of sound has not been added here yet and should be done below
                        durInd=durInd + Duration;
                    end

                    % calculate duration in sec and save the begining and end
                    % indices of the section if longer than 50ms (30ms of
                    % sound)
                    durVal = (durInd)/samprate;
                    if durInd>ceil(0.05*samprate)%section needs to be longer than 50ms (20ms of silence + 30ms of sound) 
                        nsections=nsections+1;%increment the section number for that sound
                        if endInd>length(sound_in) % the stim end by a vocalization ensure that we keep the tail of psth by extending the wavsection
                            WavIndices_temp{nsections}=[begInd length(sound_in)];
                        else
                            WavIndices_temp{nsections}=[begInd endInd];
                        end
                    end
                end
                laston=on(i);
            %end
            
        end
        
        %% Isolate first vocalization section within a window of size Lsection
        n1section=n1section+1; %increment the counter for first sections isolated
        PSTH_BG{n1section}=psthbg; %save the psth of the correspoonding background calculated earlier
        Trials_BG{n1section}=BG_Trials_temp; %save the spike arrival times of the background window on which is calculated the previous psth
        
        % Find how much silence there is between the first and second
        % vocalization sections or the end of the stim
        if nsections==2
            Sil = WavIndices_temp{nsections}(1) - WavIndices_temp{nsections-1}(2);
        elseif i==length(on)%section ends at the end of the stim
            Sil = 0; %no wave without sound available after the sectio
        else
            fprintf(1,'pb with the loop ln 341\n');
        end
        
        %define the begining and end of the section and making sure that...
        ...if the stim start by a vocalization there is still the 10ms of...
            ...silence before 
        if nsections==1 %une seule section découpée
            begInd= WavIndices_temp{nsections}(1);
            endInd= WavIndices_temp{nsections}(2);
        elseif nsections==2 %2 sections découpées
            begInd=WavIndices_temp{nsections-1}(1);
            endInd=WavIndices_temp{nsections-1}(2);
        else
            fprintf('problem with section cutting, the code cut %d when he should cut only the first 2\n', nsections);
        end
            
        if begInd==1 %the stim start by a vocalization
            prezeros=zeros(Duration,1);%padd with zero before for 10ms because was not done earlier
            if pl>0
                fprintf('the stim start by a vocalization, index of stim begining: %d\n', begInd')
            end
            begInd=-(length(prezeros)-1);%-1 because there is the zero index!
        end
            
        durInd =  endInd - begInd;
        
        
        %depending on the extracts length either enlarge the sound extract or cut
        if durInd<=Lsection
            % Extend wav file as much as possible around the vocalization
            % to include surrounding silence periods
            Needzeros = Lsection-durInd;
            if Needzeros<=Sil
                endInd_spike=endInd + Needzeros;
                begInd_spike=begInd;
            elseif Sil==0
                endInd_spike=endInd + ceil(Needzeros/2);
                begInd_spike=begInd - floor(Needzeros/2);
            else
                endInd_spike=endInd + Sil;
                begInd_spike=begInd - (Needzeros-Sil);
            end
            % write WavIndices with the new limits and save the section
            % of sound completed with needed zeros before the stim
            if begInd_spike<=0 %window begins before the start of the stim
                WavIndices{n1section}(1)=1;
                prezeros=zeros(abs(begInd_spike),1);
            else
                WavIndices{n1section}(1)=begInd_spike;
                prezeros=[];
            end
            if Sil==0 %window stops after the end of the stim
                WavIndices{n1section}(2)=length(sound_in);
                postzeros=zeros((endInd_spike - length(sound_in)),1);
            else
                WavIndices{n1section}(2)=endInd_spike;
                postzeros=[];
            end
             sec=[prezeros; sound_in(WavIndices{n1section}(1):WavIndices{n1section}(2)); postzeros];
            
            
            if length(sec)~=Lsection
                fprintf(1,'pb withduration of section, the sound length is %d when it should be %d\n', length(sec), Lsection);
            end
            
            SectionWave{n1section} = sec;
            SectionLength(n1section) = length(sec)/samprate*1000;%durée de l'extrait en ms ici
            Section_cat{n1section}='full';
        else
            %cut wav file and store
            begInd_spike=begInd;
            endInd_spike=begInd + Lsection;
            if begInd_spike<=0 %window begins before the start of the stim
                WavIndices{n1section}=[1 endInd_spike];
                prezeros=zeros(abs(begInd_spike),1);
                sec=[prezeros; sound_in(1:endInd_spike)];
            else
                WavIndices{n1section}=[begInd_spike endInd_spike];
                sec=sound_in(begInd_spike : endInd_spike);
            end
            SectionWave{n1section} = sec;
            SectionLength(n1section) = length(sec)/samprate*1000;%durée de l'extrait en ms ici
            Section_cat{n1section}='cut';
        end
                    
        %% calculate and store spectro
        % Parameters for the Spectrogram
        nstd = 6;
        fband = 50;
        twindow = 1000*nstd/(fband*2.0*pi);           % Window length in ms - 6 times the standard dev of the gaussian window
        winLength = fix(twindow*samprate/1000.0);  % Window length in number of points
        winLength = fix(winLength/2)*2;            % Enforce even window length
        increment = fix(0.001*samprate);           % Sampling rate of spectrogram in number of points - set at 1 kHz
        %calculate spectro
        [s, to, fo, pg] = GaussianSpectrum(SectionWave{n1section}, increment, winLength, samprate);
        %reshape and store spectro
        D=length(to);
        F=length(fo);
        VectorS=reshape(abs(s),1,F*D);
        Spectroto= to;
        Spectrofo= fo;
        Spectro{n1section}=VectorS;

        %% Isolate spikes that relate to the section and...
        ...calculate average (psth) for this section.
        [Section_Trials,psth,stim_mean_rate, spike_std_rate, spike_zscore, pvalue, tvalue] = spikeTimes_psth_cal(begInd_spike, endInd_spike, samprate,response,spike_rate_bg, Win);
        Trials{n1section}=Section_Trials;
        PSTH{n1section}=psth;

        %% Plot sound pressure waveform, spectrogram and psth isolated
        if pl>0
            figure(3)
            
            subplot(3,1,1);
            DBNOISE = 40;
            logB = 20*log10(abs(s));
            maxB = max(max(logB));
            minB = maxB-DBNOISE;            
            imagesc(to,fo,logB);          % to is in seconds
            axis xy;
            caxis('manual');
            caxis([minB maxB]); 
            v_axis(1) = 0;
            v_axis(2) = SectionLength(n1section)/1000;
            v_axis(3)= 0; 
            v_axis(4)= 12000;
            axis(v_axis);                                
            ylabel('Frequency kHz');
            cmap = spec_cmap();
            colormap(cmap);

            subplot(3,1,2);
            cla;
            wind1 = hanning(31)/sum(hanning(31));   % 31 ms smoothing
            smpsth = conv(psth,wind1);
            plot(((1:length(psth))/1000),smpsth(16:length(smpsth)-15)*1000);
            axis([v_axis(1) v_axis(2) 0 1000*max(smpsth)]);
            ylabel('Rate (spikes/s)');
            xlabel('Time (s)');

            subplot(3,1,3)
            cla;
            wave=SectionWave{n1section};
            plot(wave)
            pause
        end
        
        %% Store other infos on section
        VocType{n1section}=voctype{isound};
        TDT_wavfiles{n1section}=response.tdt_wavfile;
        Original_wavfiles{n1section}=response.original_wavfile;
        if strcmp(response.stim_type,'call')
            ESex{n1section}=response.stim_source_sex;
            Eage{n1section}=response.callerAge;
            Erelated{n1section}=response.stim_source;
            if ~isempty(OldWav)
                VocBank_Wavfiles{n1section}=OldWav{VocBank_idx, 3};
            end
        else
            ESex{n1section}='NaN';
            Eage{n1section}='NaN';
            Erelated{n1section}='NaN';
            if ~isempty(OldWav)
                VocBank_Wavfiles{n1section}=OldWav{VocBank_idx, 3};
            end
        end

        % Stuff values
        %iintensity2 = [iintensity2 intVal];
        histDur2 = [histDur2 durVal];
        idur2 = [idur2 durVal];
    end
    %intensity{isound}=iintensity2;
    dur{isound}=idur2;
    hist(histDur2,60)

    Res.VocType=VocType(1:n1section); % this is the type of vocalization (e.g. distance call DC, Nest call Ne, Aggressive call Ag...)
    Res.Section_cat=Section_cat(1:n1section); % identify whether the section contain an entire vocalization (full) a portion of a vocalization (cut) just the end of a vocalization (cut_end)
    Res.VocBank_wavfiles=VocBank_Wavfiles(1:n1section); % name of the wav file of the vocalization bank to which this section responded
    Res.Original_wavfiles=Original_wavfiles(1:n1section); % The real stim is a combination of 1 or 3 calls or 2.5s song. This is the original name of the wav file JEE constructed with the vocalization from the vocalization bank.
    Res.TDT_wavfiles=TDT_wavfiles(1:n1section); % The name of same previous wav stim given by TDT (stim1, stim2.... stim136...)
    Res.WavIndices = WavIndices(1:n1section); % begining and end indices of each section within the TDT_wavfiles, just to be able to retrieve the wavform if needed
    %Res.SectionWave=SectionWave(1:nsections);
    Res.SectionLength = SectionLength(1:n1section); % duration of each section in ms
    Res.ESex=ESex(1:n1section); % Sex of the emitter of the vocalization
    Res.Eage=Eage(1:n1section);  % Age of the emitter of the vocalization
    Res.Erelated=Erelated(1:n1section); % Relation of the emitter to the subject (familiar, unfamiliar, self)
    Res.Trials=Trials(1:n1section); % Contains the spike arrival times in ms from the begining of the section and not in ms from the begining of the stim as in h5 files!!!
    Res.PSTH=PSTH(1:n1section);
    Res.Spectro=Spectro(1:n1section);
    Res.Spectroto=Spectroto;
    Res.Spectrofo=Spectrofo;
    Res.VocDuration=dur;
    Res.PSTH_BG=PSTH_BG(1:n1section);
    Res.Trials_BG=Trials_BG(1:n1section);

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
        filename=fullfile('/auto','k6','julie','matfile',Res.subject,['FirstVoc_' Res.Site '.mat']);
    end
    %if sum(ttest==1)>0    
        save(filename, '-struct', 'Res');
        fprintf('saved data under: %s\n', filename);
    %else
    %    fprintf('none significant site\n');
    %end
    clear duration Res VocType TDT_wavfiles Cut_orders Original_wavfiles WavIndices SectionWave ESex Eage Erelated Trials PSTH MeanRate StdRate Spectro Section_cat Section_zscore SectionLength Section_tvalue Section_pvalue sections_good_zscores
else
    fprintf(1, 'Warning: could not find stimType %s in h5 file %s\n', stimType, h5Path);
    fprintf(1, '\tAvailable options are:\n');
    for ic=1:nclasses
        fprintf(1,'%s\n', unit.classes{ic});
    end
    clear duration unit pthreshold
end
end
    