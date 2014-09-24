function [filename]=Matfile_construct_JEE(h5Path, OldWav, pl)
%% This function cut tdt stim to isolate vocalizations as much as possible...
...and create 200ms windows of analysis.

duration=0.02; % duration in ms of min silence period
if nargin<3
    pl=0; %set to 0 for no graph; 1 if you want to see final results per stim; 2 if you to see both steps graph
end
%if nargin<2
    pthreshold=0.05;
%end

%%  I choose 200ms as a size to cut my stims
Win=0.2;

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
    Cut_orders=zeros(nfiles*coeff,1);
    Voc_orders=zeros(nfiles*coeff,1);
    Original_wavfiles=cell(nfiles*coeff,1);
    SectionWave=cell(nfiles*coeff,1);
    WavIndices=cell(nfiles*coeff,1);
    SectionLength = zeros(nfiles*coeff,1);
    ESex=cell(nfiles*coeff,1);
    Eage=cell(nfiles*coeff,1);
    Erelated=cell(nfiles*coeff,1);
    Trials=cell(nfiles*coeff,1);
    PSTH=cell(nfiles*coeff,1);
    MeanRate=cell(nfiles*coeff,1);
    MeanRateVoc=cell(nfiles*coeff,1);
    StdRate=cell(nfiles*coeff,1);
    StdRateVoc=cell(nfiles*coeff,1);
    Section_zscore=cell(nfiles*coeff,1);
    Section_tvalue=cell(nfiles*coeff,1);
    Section_pvalue=cell(nfiles*coeff,1);
    Spectro=cell(nfiles*coeff,1);
    Section_cat=cell(nfiles*coeff,1);
    histDur2=[];
    nsections = 0;
    Ag = zeros(1,10);
    Ag(1)=1;
    Be = zeros(1,10);
    Be(2)=1;
    DC = zeros(1,10);
    DC(3)=1;
    Di = zeros(1,10);
    Di(4)=1;
    LT = zeros(1,10);
    LT(5)=1;
    Ne = zeros(1,10);
    Ne(6)=1;
    Te = zeros(1,10);
    Te(7)=1;
    Th = zeros(1,10);
    Th(8)=1;
    mlnoise = zeros(1,10);
    mlnoise(9)=1;
    song = zeros(1,10);
    song(10)=1;
    sections_good_zscores=zeros(1,10);
    total_sections=zeros(1,10);


    for isound = 1:nfiles
        response=responses{isound};
        stim_name=response.tdt_wavfile;
        stim_number=str2double(response.number);
        if ~isempty(OldWav)
            VocBank_idx= find(cell2mat(OldWav(:,1))==stim_number);
        end
    
       %Set up the counting matrices depending on the type of stim
        if strcmp(voctype{isound}, 'Ag')
            MC=Ag;
        elseif strcmp(voctype{isound}, 'Be')
            MC=Be;
        elseif strcmp(voctype{isound}, 'DC')
            MC=DC;
        elseif strcmp(voctype{isound}, 'Di')
            MC=Di;
        elseif strcmp(voctype{isound}, 'LT')
            MC=LT;
        elseif strcmp(voctype{isound}, 'Ne')
            MC=Ne;
        elseif strcmp(voctype{isound}, 'Te')
            MC=Te;
        elseif strcmp(voctype{isound}, 'Th')
            MC=Th;
        elseif strcmp(voctype{isound}, 'mlnoise')
            MC=mlnoise;
        elseif strcmp(voctype{isound}, 'song')
            MC=song;
        else
            MC=zeros(1,10);
            fprintf('Stim %d of category %s not taken into account for the test of significant unit', isound, voctype{isound});
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
        iSound2 = zeros(0);
        iintensity=intensity{isound};
        Thresh=0.2*mean(iintensity); %silence is defined as 20% of average intensity calculated above
        yy = 1;
        Duration = ceil(0.06*samprate);% calls have to be spaced by at least 60ms to be separated. 60ms tail is added after the end sensus stricto of each sound
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
                on((Sil(1,cc)+Duration):1:(Sil(2,cc)-ceil(0.01*samprate)))=0;% indicate zones of silences (0) and zones around calls that should be kept (60ms after the end and 10ms before the begining)
            end
            if Sil(2,c)==length(on)%if the wavfile ends by a silence
                on((Sil(2,c)-ceil(0.01*samprate)):1:end)=0;% then make sure the silence was not supress at the end by the previous statement
            end
        end
        if on(end)==1
            on= [on;ones(Duration,1)];%in case the stim end by a vocalization, add some space after (60ms) to keep the tail of the spike rate
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
        for it=1:ntrials
            trial = response.trials{it};
            bgdur(it)=trial.pre_time;
            spike_time_trials = trial.spikeTimes;
            spike_time1=(spike_time_trials<0);
            spike_time2=(spike_time_trials>=(-bgdur(it)));
            SpikesIndices=find(spike_time1 .* spike_time2);
            spike_count_bg(it)=length(SpikesIndices);
            
        end
        %bg_mean_rate = mean(spike_count_bg)./(bgdur);
        %bg_std_rate = std(spike_count_bg)./(bgdur);
        spike_rate_bg = spike_count_bg./bgdur;
    
    
        %% Find beggining and end of each vocalization-section
        laston = 0;
        voc_ord=0;
        for i=1:length(on)
            if (laston == 0) && (on(i)~=0)
                begInd = i;
                voc_ord=voc_ord + 1; %increment the vocalization order index
            elseif (laston~=0) && (on(i) == 0) || (laston~=0) && (i==length(on))
                endInd = i;
                durInd = endInd - begInd;
                    
                % calculate duration in sec
                durVal = (durInd)/samprate;
                
                %initialize the counter of the section order within each
                %vocalization
                ord=0;
          
                %% depending on section length enlarged or regularly cut to obtain sections of 200ms each
                if durInd<=Lsection && durInd>ceil(0.03*samprate)%section shorter than or equal to 200ms
                    nsections=nsections+1;
                    ord=ord+1;
                    total_sections=total_sections+MC;
                    % Complete wav file with zeros and store 
                    if endInd>length(sound_in)
                        sec=sound_in(begInd:end);
                        WavIndices{nsections}=[begInd length(sound_in)];
                    else
                        sec=sound_in(begInd:endInd);
                        WavIndices{nsections}=[begInd endInd];
                    end
                    blanck=zeros((Lsection-length(sec)),1);
                    SectionWave{nsections} = [sec; blanck];
                    SectionLength(nsections) = Lsection/samprate*1000;
                    
                    %% calculate and store spectro
                    % Parameters for the Spectrogram
                    nstd = 6;
                    fband = 50;
                    twindow = 1000*nstd/(fband*2.0*pi);           % Window length in ms - 6 times the standard dev of the gaussian window
                    winLength = fix(twindow*samprate/1000.0);  % Window length in number of points
                    winLength = fix(winLength/2)*2;            % Enforce even window length
                    increment = fix(0.001*samprate);           % Sampling rate of spectrogram in number of points - set at 1 kHz
                    %calculate spectro
                    [s, to, fo, pg] = GaussianSpectrum(SectionWave{nsections}, increment, winLength, samprate);
                    %reshape and store spectro
                    D=length(to);
                    F=length(fo);
                    VectorS=reshape(abs(s),1,F*D);
                    Spectroto= to;
                    Spectrofo= fo;
                    Spectro{nsections}=VectorS;
                    
                    %% Isolate spikes that relate to the section and...
                    ...calculate average (psth) for this section. The psth isolated...
                        ...can be shorter than 200ms its lentgh depends on...
                        ...the length of the section and then padded with...
                        ... zeros to reach the length of the Window Win
                    [Section_Trials,psth,stim_mean_rate, spike_std_rate, spike_zscore, pvalue, tvalue] = spikeTimes_psth_cal(begInd, endInd, samprate,response,spike_rate_bg, Win);
                    Trials{nsections}=Section_Trials;
                    PSTH{nsections}=psth;
                    MeanRate{nsections}=stim_mean_rate;
                    MeanRateVoc{nsections}=stim_mean_rate;% here the spike rate is the same for the section and for the whole vocalization
                    StdRate{nsections}=spike_std_rate;
                    StdRateVoc{nsections}=spike_std_rate;
                    Section_zscore{nsections}=spike_zscore;
                    Section_tvalue{nsections}=tvalue;
                    Section_pvalue{nsections}=pvalue;
                    if pvalue<pthreshold
                        sections_good_zscores=sections_good_zscores+MC;
                    end
                    
                    %% Store other infos on section
                    VocType{nsections}=voctype{isound};
                    TDT_wavfiles{nsections}=response.tdt_wavfile;
                    Cut_orders(nsections)=ord;
                    Voc_orders(nsections)=voc_ord;
                    Original_wavfiles{nsections}=response.original_wavfile;
                    if strcmp(response.stim_type,'call')
                        ESex{nsections}=response.stim_source_sex;
                        Eage{nsections}=response.callerAge;
                        Erelated{nsections}=response.stim_source;
                        if ~isempty(OldWav)
                            if ~strcmp(response.callid, 'Be')
                                VocBank_Wavfiles{nsections}=OldWav{VocBank_idx, voc_ord+2};
                            else
                                VocBank_Wavfiles{nsections}=OldWav{VocBank_idx, 3};
                            end
                        end
                    else
                        ESex{nsections}='NaN';
                        Eage{nsections}='NaN';
                        Erelated{nsections}='NaN';
                        if ~isempty(OldWav)
                            VocBank_Wavfiles{nsections}=OldWav{VocBank_idx, 3};
                        end
                    end
                    Section_cat{nsections}='full';
                
                elseif durInd>Lsection %section longer than 200ms
                    % Calculate the mean rate and tsd rate of that long
                    % section
                    [stim_mean_rateVoc, spike_std_rateVoc, stim_rateVoc] = meanrate_cal(begInd, endInd, samprate,response);
                    
                    nsec=floor(durInd/Lsection); %find the number of 200ms sections we can have from this long section
                    for ns=1:nsec
                        ord=ord+1;
                        nsections=nsections+1;
                        total_sections=total_sections+MC;
                        
                        % Cut and store wav file
                        if (begInd+ns*Lsection)>length(sound_in)
                            local_begInd = (begInd+(ns-1)*Lsection);
                            local_endInd = length(sound_in);
                            sec=sound_in(local_begInd:end,1);
                            WavIndices{nsections}=[local_begInd local_endInd];
                            blanck=zeros((Lsection-length(sec)),1);
                            SectionWave{nsections} = [sec; blanck];
                            SectionLength(nsections) = Lsection/samprate*1000;
                        else
                            local_begInd = (begInd+(ns-1)*Lsection);
                            local_endInd = (begInd+ns*Lsection);
                            sec=sound_in(local_begInd:local_endInd,1);
                            WavIndices{nsections}=[local_begInd local_endInd];
                            SectionWave{nsections} = sec;
                            SectionLength(nsections) = length(sec)/samprate*1000;
                        end
                    
                        % calculate and store spectro
                        % Parameters for the Spectrogram
                        nstd = 6;
                        fband = 50;
                        twindow = 1000*nstd/(fband*2.0*pi);           % Window length in ms - 6 times the standard dev of the gaussian window
                        winLength = fix(twindow*samprate/1000.0);  % Window length in number of points
                        winLength = fix(winLength/2)*2;            % Enforce even window length
                        increment = fix(0.001*samprate);           % Sampling rate of spectrogram in number of points - set at 1 kHz
                        %calculate spectro
                        [s, to, fo, pg] = GaussianSpectrum(SectionWave{nsections}, increment, winLength, samprate);
                        %reshape and store spectro
                        D=length(to);
                        F=length(fo);
                        VectorS=reshape(abs(s),1,F*D);
                        Spectroto= to;
                        Spectrofo= fo;
                        Spectro{nsections}=VectorS;
                    
                        % Isolate spikes that relate to the section and calculate average (psth) for this section
                        [Section_Trials,psth,stim_mean_rate, spike_std_rate, spike_zscore, pvalue, tvalue] = spikeTimes_psth_cal(local_begInd, local_endInd, samprate,response,spike_rate_bg, Win);
                        Trials{nsections}=Section_Trials;
                        PSTH{nsections}=psth;
                        MeanRate{nsections}=stim_mean_rate;%Mean rate of the section (200ms) = slice of vocalization
                        MeanRateVoc{nsections}=stim_mean_rateVoc;%Mean rate of the isolated vocalization before slicing
                        StdRate{nsections}=spike_std_rate;%Std rate of the section (200ms) = slice of vocalization
                        StdRateVoc{nsections}=spike_std_rateVoc;%Std rate of the isolated vocalization before slicing
                        Section_zscore{nsections}=spike_zscore;
                        Section_tvalue{nsections}=tvalue;
                        Section_pvalue{nsections}=pvalue;
                        if pvalue<pthreshold
                            sections_good_zscores=sections_good_zscores+MC;
                        end
                
                        % Store other infos on section
                        VocType{nsections}=voctype{isound};
                        TDT_wavfiles{nsections}=response.tdt_wavfile;
                        Cut_orders(nsections)=ord;
                        Voc_orders(nsections)=voc_ord;
                        Original_wavfiles{nsections}=response.original_wavfile;
                        if strcmp(response.stim_type,'call')
                            ESex{nsections}=response.stim_source_sex;
                            Eage{nsections}=response.callerAge;
                            Erelated{nsections}=response.stim_source;
                            if ~isempty(OldWav)
                                if ~strcmp(response.callid, 'Be')
                                    VocBank_Wavfiles{nsections}=OldWav{VocBank_idx, voc_ord+2};
                                else
                                    VocBank_Wavfiles{nsections}=OldWav{VocBank_idx, 3};
                                end
                            end
                        else
                            ESex{nsections}='NaN';
                            Eage{nsections}='NaN';
                            Erelated{nsections}='NaN';
                            if ~isempty(OldWav)
                                VocBank_Wavfiles{nsections}=OldWav{VocBank_idx, 3};
                        
                            end
                        end
                        Section_cat{nsections}='cut';
                    end
                    % Treat the end of the wavsection if its longer than
                    % 100ms (40ms call + 60ms of forced silence) otherwise disgard
                    if (durInd-nsec*Lsection)>0.1*samprate
                        % Cut and store wav file
                        ord=ord+1;
                        nsections=nsections+1;
                        total_sections=total_sections+MC;
                        if endInd>length(sound_in)
                            local_begInd=begInd+nsec*Lsection;
                            local_endInd=length(sound_in);
                            sec=sound_in(local_begInd:end,1);
                            WavIndices{nsections} = [local_begInd local_endInd];
                        else
                            local_begInd=begInd+nsec*Lsection;
                            local_endInd=endInd;
                            sec=sound_in((begInd+nsec*Lsection):endInd,1);
                            WavIndices{nsections} = [(begInd+nsec*Lsection) endInd];
                        end
                        blanck=zeros((Lsection-length(sec)),1);
                        SectionWave{nsections} = [sec ; blanck];  
                        SectionLength(nsections) = Lsection/samprate*1000;
                    
                        % calculate and store spectro
                        % Parameters for the Spectrogram
                        nstd = 6;
                        fband = 50;
                        twindow = 1000*nstd/(fband*2.0*pi);           % Window length in ms - 6 times the standard dev of the gaussian window
                        winLength = fix(twindow*samprate/1000.0);  % Window length in number of points
                        winLength = fix(winLength/2)*2;            % Enforce even window length
                        increment = fix(0.001*samprate);           % Sampling rate of spectrogram in number of points - set at 1 kHz
                        %calculate spectro
                        [s, to, fo, pg] = GaussianSpectrum(SectionWave{nsections}, increment, winLength, samprate);
                        %reshape and store spectro
                        D=length(to);
                        F=length(fo);
                        VectorS=reshape(abs(s),1,F*D);
                        Spectro{nsections}=VectorS;
                    
                        % Isolate spikes that relate to the section and calculate average (psth) for this section
                        [Section_Trials,psth,stim_mean_rate, spike_std_rate, spike_zscore, pvalue, tvalue] = spikeTimes_psth_cal(local_begInd, local_endInd, samprate,response,spike_rate_bg, Win);
                        Trials{nsections}=Section_Trials;
                        PSTH{nsections}=psth;
                        MeanRate{nsections}=stim_mean_rate;%Mean rate of the section (200ms) = slice of vocalization
                        MeanRateVoc{nsections}=stim_mean_rateVoc;%Mean rate of the isolated vocalization before slicing
                        StdRate{nsections}=spike_std_rate;%Std rate of the section (200ms) = slice of vocalization
                        StdRateVoc{nsections}=spike_std_rateVoc;%Std rate of the isolated vocalization before slicing
                        Section_zscore{nsections}=spike_zscore;
                        Section_tvalue{nsections}=tvalue;
                        Section_pvalue{nsections}=pvalue;
                        if pvalue<pthreshold
                            sections_good_zscores=sections_good_zscores+MC;
                        end
                    
                        % Store other infos on section
                        VocType{nsections}=voctype{isound};
                        TDT_wavfiles{nsections}=response.tdt_wavfile;
                        Cut_orders(nsections)=ord;
                        Voc_orders(nsections)=voc_ord;
                        Original_wavfiles{nsections}=response.original_wavfile;
                        if strcmp(response.stim_type,'call')
                            ESex{nsections}=response.stim_source_sex;
                            Eage{nsections}=response.callerAge;
                            Erelated{nsections}=response.stim_source;
                            if ~isempty(OldWav)
                                if ~strcmp(response.callid, 'Be')
                                    VocBank_Wavfiles{nsections}=OldWav{VocBank_idx, voc_ord+2};
                                else
                                    VocBank_Wavfiles{nsections}=OldWav{VocBank_idx, 3};
                                end
                            end
                        else
                            ESex{nsections}='NaN';
                            Eage{nsections}='NaN';
                            Erelated{nsections}='NaN';
                            if ~isempty(OldWav)
                                VocBank_Wavfiles{nsections}=OldWav{VocBank_idx, 3};
                            end
                        end
                        Section_cat{nsections}='cut_end';
                    end
                end
                % Stuff values
                %iintensity2 = [iintensity2 intVal];
                histDur2 = [histDur2 durVal];
                idur2 = [idur2 durVal];
                iSound2 = [iSound2 [begInd ; endInd]];
                % fprintf(1, 'Beg= %.2f End= %.2f Int= %.2f Dur= %.2f Freq= %.0f RD= %.2f\n', begTime, endTime, 20*log10(intVal), durVal, fVal, rateDiffVal);
            end
            laston = on(i);
        end
        %intensity{isound}=iintensity2;
        dur{isound}=idur2;
        Sound{isound}=iSound2;
    end
    hist(histDur2,60)

    Res.VocType=VocType(1:nsections); % this is the type of vocalization (e.g. distance call DC, Nest call Ne, Aggressive call Ag...)
    Res.Voc_orders=Voc_orders(1:nsections); % identify the vocalization from which the section was isolated. Because they are often 3 calls per tdt stim or a lot of syllables, each section can be cut from a different vocalization within the same tdt stim
    Res.Cut_orders=Cut_orders(1:nsections); % identify the section position within the same vocalization within the same tdt stim
    Res.Section_cat=Section_cat(1:nsections); % identify whether the section contain an entire vocalization (full) a portion of a vocalization (cut) just the end of a vocalization (cut_end)
    Res.VocBank_wavfiles=VocBank_Wavfiles(1:nsections); % name of the wav file of the vocalization bank to which this section responded
    Res.Original_wavfiles=Original_wavfiles(1:nsections); % The real stim is a combination of 1 or 3 calls or 2.5s song. This is the original name of the wav file JEE constructed with the vocalization from the vocalization bank.
    Res.TDT_wavfiles=TDT_wavfiles(1:nsections); % The name of same previous wav stim given by TDT (stim1, stim2.... stim136...)
    Res.WavIndices = WavIndices(1:nsections); % begining and end indices of each section within the TDT_wavfiles, just to be able to retrieve the wavform if needed
    %Res.SectionWave=SectionWave(1:nsections);
    Res.SectionLength = SectionLength(1:nsections); % duration of each section in ms
    Res.ESex=ESex(1:nsections); % Sex of the emitter of the vocalization
    Res.Eage=Eage(1:nsections);  % Age of the emitter of the vocalization
    Res.Erelated=Erelated(1:nsections); % Relation of the emitter to the subject (familiar, unfamiliar, self)
    Res.Trials=Trials(1:nsections); % Contains the spike arrival times in ms from the begining of the section and not in ms from the begining of the stim as in h5 files!!!
    Res.PSTH=PSTH(1:nsections);
    Res.MeanRate=MeanRate(1:nsections);
    Res.StdRate=StdRate(1:nsections);
    Res.MeanRateVoc=MeanRateVoc(1:nsections);
    Res.StdRateVoc=StdRate(1:nsections);
    Res.Spectro=Spectro(1:nsections);
    Res.Spectroto=Spectroto;
    Res.Spectrofo=Spectrofo;
    Res.Section_zscore=Section_zscore(1:nsections);
    Res.Section_tvalue=Section_tvalue(1:nsections);
    Res.Section_pvalue=Section_pvalue(1:nsections);
    Res.sections_good_zscores=sections_good_zscores;
    Res.total_sections=total_sections;


    %save all the results for this cell in a mat file only if the cell always
    %respond (pvalue<pthreshold) to at least one category of stim
    ttest=sections_good_zscores./total_sections;

    if ismac()
            [status username] = system('who am i');
            if strcmp(strtok(username), 'frederictheunissen')
                if strncmp('/auto/fdata/solveig',stim_name, 19)
                elseif strncmp('/auto/fdata/julie',stim_name, 17)
                    filename = fullfile('/Users','frederictheunissen','Documents','Data','Julie','matfile',Res.subject,[Res.Site '.mat']);
                end
            elseif strcmp(strtok(username), 'elie')
                filename = fullfile('/Users','elie','Documents','MATLAB','data','matfile',Res.subject,['ZS_' Res.Site '.mat']);
            end
    else
        filename=fullfile('/auto','k6','julie','matfile',Res.subject,['ZS_' Res.Site '.mat']);
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
    