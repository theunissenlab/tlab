function [Res]=VocSectioning(Folder, pl)
%% This function cut sound extract to isolate vocalizations as much as possible, and create 200ms windows of analysis.
...it's running in twe steps. The first step calculate the mean intensity of quickly localized sound objects.
    ...the second step refines each section of sound by substracting all the silence periods (all the periods
    ...where the sound intensity is below 20% of the mean sound intensity calculated at the first step).
    ...WARNING: this code adds 60ms of silence after each

%Set the parameters for sound filtering 
songfilt_args.f_low = 250;
songfilt_args.f_high = 12000;
songfilt_args.db_att = 0;
duration=0.02; % duration in ms of min silence period

if nargin<2
    pl=0; %set to 0 for no graph; 1 if you want to see final results per stim; 2 if you to see both steps graph
end

if nargin<1
    %% Find the data base on the cluster on a local mac machine.
    if ismac()
        [status username] = system('who am i');
        if strcmp(strtok(username), 'frederictheunissen')
            DataBasePath='/Users/frederictheunissen/Documents/Data/Julie';
        elseif strcmp(strtok(username), 'elie')
            DataBasePath='/Users/elie/Documents/ManipBerkeley/Recordings/group recordings';
        end
    end
else
    DataBasePath=Folder;
SoundFiles=dir(fullfile(DataBasePath, '*.wav'));

%% Cutting procedure based on zero values to calculate intensity values of each cut
    % This is the number of sound files
    nfiles = size(SoundFiles,1);
    
    % set up output cell arrays
    intensity = cell(nfiles,1);
    dur = cell(nfiles,1);
    voctype = cell(nfiles,1);
    Sound = cell(nfiles,1);
    histDur=[];

    %Loop through files, detect sound periods
    for isound = 1:nfiles
        stim_name=SoundFiles(isound).name;
        [sound_in, samprate] = wavread(fullfile(DataBasePath, stim_name));
        
        % Read the callid.
        [SoundPath, SoundName, ext]=fileparts(stim_name);
        vocid=SoundName(19:20);
        if size(sound_in,2)==2
            sound_in=(sound_in(:,1) + sound_in(:,2))/2;
        end
        sound_in = songfilt_call(sound_in,samprate,songfilt_args.f_low,songfilt_args.f_high,...
                 songfilt_args.db_att, vocid);
        iintensity=[];
        idur=[];
        iSound = zeros(0);
        
    
       
        
    
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
                plot_sound_selection(vocid, sound_in, samprate, on);
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
    %hist(histDur,60)



    %%  I choose 200ms as a size to cut my stim
    Lsection=ceil(0.2*samprate);



    %% Now run a second cutting system based on the intensity of the signal to get nicer cuts and enlarged the saved window around the vocalization
    
    % prepare cell arrays
    coeff=8; %estimate of 8 sections per stim on average
    VocType=cell(nfiles*coeff,1);
    Cut_orders=zeros(nfiles*coeff,1);
    Voc_orders=zeros(nfiles*coeff,1);
    Original_wavfiles=cell(nfiles*coeff,1);
    SectionWave=cell(nfiles*coeff,1);
    WavIndices=cell(nfiles*coeff,1);
    SectionLength = zeros(nfiles*coeff,1);
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
    Wh = zeros(1,10);
    Wh(9)=1;
    So = zeros(1,10);
    So(10)=1;
    total_sections=zeros(1,10);


    for isound = 1:nfiles
        ord=0;
        stim_name=SoundFiles(isound).name;
        [sound_in, samprate] = wavread(fullfile(DataBasePath, stim_name)); 
        if size(sound_in,2)==2
            sound_in=(sound_in(:,1) + sound_in(:,2))/2;
        end
        sound_in = songfilt_call(sound_in,samprate,songfilt_args.f_low,songfilt_args.f_high,...
                 songfilt_args.db_att, voctype{isound});
    
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
        elseif strcmp(voctype{isound}, 'So')
            MC=So;
        elseif strcmp(voctype{isound}, 'Wh')
            MC=Wh;
        else
            MC=zeros(1,10);
            fprintf('Stim %d of category %s not taken into account for the test of significant unit', isound, voctype{isound});
        end
        
        %Find Silences longer than 60ms between vocalizations
        idur2=[];
        iSound2 = zeros(0);
        iintensity=intensity{isound};
        Thresh=0.2*mean(iintensity); %silence is defined as 40% of mean intensity calculated above
        yy = 1;
        Duration = ceil(0.06*samprate);% vocalizations have to be spaced by at least 60ms to be separated. 60ms tail is added after the end sensus stricto of each sound
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
                on((Sil(1,cc)+ceil(0.01*samprate)):1:(Sil(2,cc)-ceil(0.01*samprate)))=0;% indicate zones of silences (0) and zones around calls that should be kept (60ms after the end and 10ms before the begining)
            end
            if Sil(2,c)==length(on)%if the wavfile ends by a silence
                on((Sil(2,c)-ceil(0.01*samprate)):1:end)=0;% then make sure the silence was not supress at the end by the previous statement
            end
        end
        if on(end)==1
            on= [on;ones(ceil(0.01*samprate),1)];%in case the stim end by a vocalization, add some space after (10ms)
        end
        on=on*0.8;
        if pl>0        
            plot_sound_selection(voctype{isound}, sound_in, samprate, on);
            pause;
        end
    
        %% Find beggining and end of each vocalization-section
        laston = 0;
        voc_ord=0;
        for i=1:length(on)
            if (laston == 0) && (on(i)~=0)
                begInd = i;
                voc_ord=voc_ord+1;
            elseif (laston~=0) && (on(i) == 0) || (laston~=0) && (i==length(on))
                endInd = i;
                durInd = endInd - begInd;
                    
                % calculate duration in sec
                durVal = (durInd)/samprate;
          
                %% depending on section length enlarged or regularly cut to obtain sections of 200ms each
                if durInd<Lsection && durInd>ceil(0.03*samprate)%section shorter than 200ms but longer than 30ms
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
                    
                    %% Store info on section
                    VocType{nsections}=voctype{isound};
                    Cut_orders(nsections)=ord;
                    Voc_orders(nsections)=voc_ord;
                    Original_wavfiles{nsections}=stim_name;
                    Section_cat{nsections}='full';
                
                elseif durInd>Lsection %section longer than 200ms
                    nsec=floor(durInd/Lsection); %find the number of 200ms sections we can have from this long section
                    for ns=1:nsec
                        ord=ord+1;
                        nsections=nsections+1;
                        total_sections=total_sections+MC;
                        
                        %% Cut and store wav file
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

                        %% Store info on section
                        VocType{nsections}=voctype{isound};
                        Cut_orders(nsections)=ord;
                        Voc_orders(nsections)=voc_ord;
                        Original_wavfiles{nsections}=stim_name;
                        if ns==1
                            Section_cat{nsections}='cut-first';
                        else
                            Section_cat{nsections}='cut-middle';
                        end
                    end
                    % Treat the end of the wavsection if its longer than
                    % 40ms call otherwise disgard
                    if (durInd-nsec*Lsection)>0.04*samprate
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
                            WavIndices{nsections} = [local_begInd local_endInd];
                        end
                        blanck=zeros((Lsection-length(sec)),1);
                        SectionWave{nsections} = [sec ; blanck];  
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

                        % Store info on section
                        VocType{nsections}=voctype{isound};
                        Cut_orders(nsections)=ord;
                        Voc_orders(nsections)=voc_ord;
                        Original_wavfiles{nsections}=stim_name;
                        Section_cat{nsections}='cut_end';
                    end
                end
            end
            laston = on(i);
        end
        %intensity{isound}=iintensity2;
        dur{isound}=idur2;
        Sound{isound}=iSound2;
    end
    %hist(histDur2,60)
    Res=struct;
    Res.WaveName=Original_wavfiles(1:nsections);
    Res.VocalizationType=VocType(1:nsections);
    Res.VocalizationOrder=Voc_orders(1:nsections);
    Res.Section_Order=Cut_orders(1:nsections);
    Res.Section_Category=Section_cat(1:nsections);
    %Res.Section_Wave=SectionWave{1:nsections};
    Res.Section_WaveIndices=WavIndices(1:nsections);
    Res.Section_Spectro=Spectro(1:nsections);
    Res.Spectro_to=Spectroto;
    Res.Spectro_fo=Spectrofo;
    
    
    
    clear duration VocType Cut_orders Original_wavfiles WavIndices SectionWave Spectro Section_cat SectionLength

end
    