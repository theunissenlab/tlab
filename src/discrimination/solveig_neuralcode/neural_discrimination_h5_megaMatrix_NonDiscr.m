function dataMega = neural_discrimination_h5_megaMatrix_NonDiscr(h5Path, Birdname)
%[nfiles, stimrate, stdstimrate, backrate, stdbackrate, avgzscore, stdzscore, avgdprime, stddprime, percorrect, avgrank, info, noiseinfo, rateinfo, info1, info2, info2Interval]
% Calculates 3 measures of neural discrimination, the within dprime, the
% percent correct of an ideal oberver based on VR method and the
% information based on Gamma model of spiking neuron
% read all the spike arrival times in a folder
% stimType : 'Mask1' etc. Was in the function input in previous versions of this script, not anymore in the 'Solveig' version (see l.23).

%Birdname = 'GraBla1602F'; % example...
%h5Path = 'Site1_R2601_e10_s0_ss1.h5';

% wind = 3 or a vector or different window sizes

%% Read input data from data base
unit = read_unit_h5file(h5Path, 'r'); % 'r' to read only, 'a' to add fields

%% Fix some parameters
if nargin<2
    PathsepIndices=strfind(unit.source_directory, '/');
    Birdname=unit.source_directory(PathsepIndices(end)+1:end);
end
if (strcmp(Birdname, 'GreWhi2513F') || strcmp(Birdname, 'BlaLbl1986M')) % birds not tested with Syn stimuli
    isSyn = 0;
else
    isSyn = 1;
end
isSyn = 0;
%outputFolder = '/auto/k8/fdata/solveig/matfile_NonDiscrUnits_SpikeShape';
outputFolder = '/Users/elie/Dropbox/Solveig/Neuro/J. Neuroscience/Codes & data/Best Bird and Dist -units'; %for the purpose of getting VanRossum distances between stims


%% Select the protocol (SelX, Callx, Maskx, STRFx...)
%number of different protocols run for that unit, identify if the one
%you are asking for is existing in this unit and selecting this one in responses
nclasses = length(unit.classes);
if nclasses == 1
    stimTypeCell = unit.classes;
else
    fprintf('WARNING: more than one stimType (ex. Mask1) in this unit!\n');
    return
end
stimType = stimTypeCell{1};

classId = 0;
for ic=1:nclasses
    if strcmp(stimType, unit.classes{ic})
        classId = ic;
        break;
    end
end
if classId == 0
    fprintf(1, 'Warning: could not find stimType %s in h5 file %s\n', stimType, h5Path);
    fprintf(1, '\tAvailable options are:\n');
    for ic=1:nclasses
        fprintf(1,'\t\t%s\n', unit.classes{ic});
    end
    return;
end

responses = unit.class_responses.(stimType);
%nresp = size(responses,2);

%% Construct the input/output we need for all the following code
% This is the number of sound files played
nfiles = length(responses);

%This is the vector of the length of the different stims
stim_len = zeros (nfiles,1);
stim_class = cell (nfiles,1);
birdid = cell (nfiles,1);
distance = cell (nfiles,1);
% This is the cell array of the time of spike arrival
spike_times = cell(1, nfiles);
% Calculate the average spike rate over trials for each stim
spike_rate = nan(nfiles,1);

for nf=1:nfiles
    stim_len(nf) = 1000*responses{nf}.stim_duration; %here stim_len is in ms
    stim_class{nf} = responses{nf}.stim_class;
    birdid{nf} = responses{nf}.birdid;
    distance{nf} = responses{nf}.distance;
    Trials = responses{nf}.trials;
    nt = length(Trials);
    spike_times{nf} = cell(1, nt);
    for it = 1:nt
        trial = Trials{it};
        spike_times{nf}{it} = trial.spikeTimes.*1000.0;   % spike_times are in ms in neural_discrimination and in s in hdf5
    end
    spike_rate(nf) = sum((cell2mat(spike_times{nf})>0) .* (cell2mat(spike_times{nf})<=stim_len(nf)))*1000/(length(spike_times{nf})*stim_len(nf)); % number of spikes per second
end

MaxStimLen = max(stim_len);


% This is various stuff we'll need for the final outputs
ConSyn = {'Con', 'Syn'};
dataMega = cell(1, 2);      % Final data cell (Con and Syn)

%% Set up parameters for all the calculations
% Window sizes for calculations that depend on window length
%winSize = [1 3 5 10 30 50 100];    % Window size for static Info calculation and for ideal observer
% We will run the code for 'Con' and 'Syn' calls with all winSizes.
wind = [1 3 5 10 30 50 100 250 300 MaxStimLen];
ns = length(wind);

% Reorder the distance vector from 2 m to 256 m
dist=unique(distance); % !!!!! wrong order : 128, 16, 2, 256, 64 m.
di=nan(length(dist),1);
for i=1:length(dist)
    di(i) = str2double(dist{i});           % change string to number...
end
di = sort(di);                      % reorder the distances...
StimType = cell(length(dist),1);
for i=1:length(dist)
    StimType{i} = num2str(di(i));       % put it back into a cell of strings
end
NstimType=length(StimType); % = number of distances

if (isSyn == 1)
    todo = [1 2];
elseif (isSyn == 0)
    todo = 1;
end

%% Calculations!!

if nfiles
    for qq = todo % 1: Conspecific calls, 2 : Synthetic calls (if applicable)
        
        if qq==1 % Compare bird identities depending on distance, for Conspecific calls
            fprintf('Dealing with CON stims...\n');
            NDist = NstimType; % all distances for 'Con' calls
            megaData = cell(1, 2);
            %PropDist = {'2m', '16m', '64m', '128m', '256m'};
            
            % Create the data set for each distance separately, using only Conspecific calls
            SpikeRate = cell(NDist,1);
            SpikeTrains=cell(NDist,1);
            Sections_len=cell(NDist,1);
            Stim_ID=cell(NDist, 1); % Keep track of the stim's birdid and call number, and distance
            % Stim_ID : col.1 = birdid, col. 2 = callnum, col. 3 = distance
            % the lines in this cell correspond to the lines in SpikeTrains and Sections_len
            
            for vt=1:NDist
                st=StimType(vt);
                selector=strcmp(distance, st);
                selector2=strcmp(stim_class, 'Con');
                selectorID=find(selector);
                selectorCon=find(selector2);
                selectorInd=intersect(selectorCon, selectorID);
                NselectorInd=length(selectorInd);
                SpikeTrainslocal = cell(NselectorInd,1);
                Sections_lenlocal = zeros(NselectorInd,1);
                SpikeRate_local = nan(NselectorInd,1);
                Stim_ID_local = cell(NselectorInd, 3);
                for NSI = 1:NselectorInd
                    SI=selectorInd(NSI);
                    SpikeTrainslocal{NSI} = spike_times{SI};
                    Sections_lenlocal(NSI)= stim_len(SI);
                    SpikeRate_local(NSI)=spike_rate(SI);
                    Stim_ID_local{NSI, 1} = responses{SI}.birdid;
                    Stim_ID_local{NSI, 2} = responses{SI}.callnum;
                    Stim_ID_local{NSI, 3} = responses{SI}.distance;
                end
                SpikeTrains{vt}=SpikeTrainslocal; % changer le nom ou faire une méta-structure si on fait tourner la boucle Con-Syn...
                Sections_len{vt}=Sections_lenlocal;
                Stim_ID{vt}=Stim_ID_local;
                SpikeRate{vt} = SpikeRate_local;
            end
            
            % Reorganise the cells to have all calls from the same birdid in order
            SpikeTrainsR=cell(NDist,1);
            Sections_lenR=cell(NDist,1);
            SpikeRateR=cell(NDist,1);
            Stim_IDR=cell(NDist, 1);
            birdID=unique(Stim_ID{vt}(:,1)); % même liste pour ttes distances
            NbirdID=length(birdID);
            for vt=1:NDist
                for bird=1:NbirdID
                    bi=birdID(bird);
                    select=strcmp(Stim_ID{vt}(:,1), bi);
                    selectorId=find(select);
                    NselectorId=length(selectorId);
                    incr = (bird-1)*NselectorId;
                    for ind=1:NselectorId
                        Stim_IDR{vt}{incr+ind,1} = Stim_ID{vt}{selectorId(ind),1};
                        Stim_IDR{vt}{incr+ind,2} = Stim_ID{vt}{selectorId(ind),2};
                        Stim_IDR{vt}{incr+ind,3} = Stim_ID{vt}{selectorId(ind),3};
                        SpikeTrainsR{vt}{incr+ind} = SpikeTrains{vt}{selectorId(ind)};
                        SpikeRateR{vt}(incr+ind) = SpikeRate{vt}(selectorId(ind));
                        Sections_lenR{vt}(incr+ind) = Sections_len{vt}(selectorId(ind));
                    end
                end
            end
            
            
        elseif qq==2 % Compare bird identities depending on distance, for Synthetic calls
            fprintf('Dealing with SYN stims...\n');
            NDist = NstimType-1; % All distances except 2m for 'Syn' calls
            megaData = cell(1, 2);
            %PropDist = {'16m', '64m', '128m', '256m'};
            
            SpikeTrains=cell(NDist,1);
            Sections_len=cell(NDist,1);
            Stim_ID=cell(NDist, 1); % Keep track of the stim's birdid and call number, and distance
            % Stim_ID : col.1 = birdid, col. 2 = callnum, col. 3 = distance
            
            for vt=2:NDist+1 % we search from 16 m to 256 m
                st=StimType(vt);
                selector=strcmp(distance, st);
                selector2=strcmp(stim_class, 'Syn');
                selectorID=find(selector);
                selectorSyn=find(selector2);
                selectorInd=intersect(selectorSyn, selectorID);
                NselectorInd=length(selectorInd);
                SpikeTrainslocal = cell(NselectorInd,1);
                Sections_lenlocal = zeros(NselectorInd,1);
                Stim_ID_local = cell(NselectorInd, 3);
                for NSI = 1:NselectorInd
                    SI=selectorInd(NSI);
                    SpikeTrainslocal{NSI} = spike_times{SI};
                    Sections_lenlocal(NSI)= stim_len(SI);
                    Stim_ID_local{NSI, 1} = responses{SI}.birdid;
                    Stim_ID_local{NSI, 2} = responses{SI}.callnum;
                    Stim_ID_local{NSI, 3} = responses{SI}.distance;
                end
                SpikeTrains{vt-1}=SpikeTrainslocal;
                Sections_len{vt-1}=Sections_lenlocal;
                Stim_ID{vt-1}=Stim_ID_local;
            end
            
            % Reorganise the cells to have all calls from the same birdid in order
            SpikeTrainsR=cell(NDist,1);
            Sections_lenR=cell(NDist,1);
            SpikeRateR=cell(NDist,1);
            Stim_IDR=cell(NDist, 1);
            birdID=unique(Stim_ID{vt-1}(:,1));
            NbirdID=length(birdID);
            for vt=1:NDist
                for bird=1:NbirdID
                    bi=birdID(bird);
                    select=strcmp(Stim_ID{vt}(:,1), bi);
                    selectorId=find(select);
                    NselectorId=length(selectorId);
                    incr = (bird-1)*NselectorId;
                    for ind=1:NselectorId
                        Stim_IDR{vt}{incr+ind,1} = Stim_ID{vt}{selectorId(ind),1};
                        Stim_IDR{vt}{incr+ind,2} = Stim_ID{vt}{selectorId(ind),2};
                        Stim_IDR{vt}{incr+ind,3} = Stim_ID{vt}{selectorId(ind),3};
                        SpikeTrainsR{vt}{incr+ind} = SpikeTrains{vt}{selectorId(ind)};
                        SpikeRateR{vt}(incr+ind)=SpikeRate{vt}(selectorId(ind));
                        Sections_lenR{vt}(incr+ind) = Sections_len{vt}(selectorId(ind));
                    end
                end
            end
        end
        ncalls = NselectorInd/NbirdID;
        
        %% ********Start of the calculations!***********
        
        %% for 'Con' calls, we calculate the mutual info for all window Sizes, regrouping stims per bird for all distances
        percorrectB = zeros(1,ns);
        mi_confusionB = zeros(1,ns);
        zdistTB = zeros(1,ns);              % Uses std for other for both self and other in KL Divergence calc
        pzdistTB = zeros(1,ns);
        mizdistTB = zeros(1,ns);
        mizdistncTB = zeros(1,ns);
        zdistSB = zeros(1,ns);              % Recenters "other" distrubtions to zero to minimize additive variance across songs
        pzdistSB = zeros(1,ns);             % and collapses distances across all trials of a single stimulus
        mizdistSB = zeros(1,ns);
        mizdistncSB = zeros(1,ns);
        megaMat_bird = cell(1,ns);
        ConfMat_bird = cell(1,ns);
        mi_conf_perBird = zeros(1,ns);
        percorrect_perBird = zeros(1,ns);
        VRDistance_Stims = cell(1,ns);
        VR_ResponseStims = cell(1,ns);
        mi_confVR_Bird = nan(1,ns);
        mi_confVR_Dist = nan(1,ns);
        ConfMatBirdVR =  cell(1,ns);
        ConfMatDistVR = cell(1,ns);
        
        
        % organise the megaMatrix per bird
        incr = ncalls*NDist; % nb of regrouped calls per bird
        spike_Times = cell(1, NselectorInd*NDist);
        spike_Rate = nan(1,NselectorInd*NDist);
        stim_Len = zeros(1, NselectorInd*NDist);
        verif = cell(NselectorInd*NDist, 3);
        for ii = 1:NbirdID
            Trains_temp = cell(1, incr);
            Len_temp = zeros(1, incr);
            verif_temp = cell(incr, 3);
            spike_temp = nan(1,incr);
            for jj = 1:NDist
                Trains_temp((jj-1)*ncalls+1:(jj-1)*ncalls+ncalls) = SpikeTrainsR{jj}((ii-1)*ncalls+1:(ii-1)*ncalls+ncalls);
                Len_temp((jj-1)*ncalls+1:(jj-1)*ncalls+ncalls) = Sections_lenR{jj}((ii-1)*ncalls+1:(ii-1)*ncalls+ncalls);
                verif_temp((jj-1)*ncalls+1:(jj-1)*ncalls+ncalls, :) = Stim_IDR{jj}((ii-1)*ncalls+1:(ii-1)*ncalls+ncalls, :);
                spike_temp((jj-1)*ncalls+1:(jj-1)*ncalls+ncalls) = SpikeRateR{jj}((ii-1)*ncalls+1:(ii-1)*ncalls+ncalls);
            end
            spike_Times((ii-1)*incr+1:(ii-1)*incr+(incr)) = Trains_temp;
            stim_Len((ii-1)*incr+1 : (ii-1)*incr+(incr)) = Len_temp;
            verif((ii-1)*incr+1:(ii-1)*incr+(incr), :) = verif_temp;  % pour vérification de la classif
            spike_Rate((ii-1)*incr+1 : (ii-1)*incr+(incr)) = spike_temp;
        end
        nFiles=length(spike_Times);
        
        %% Calculate the Spike Rate distance between stimuli and construct the corresponding confusion matrix
        [~, mi_confSR, ~, ~, ~, ~, ~, ~, ~, ~, ConfMatSR, ~,SRDistance_Stims, SR_ResponseStims] = info_distanceB_slvg(nFiles, spike_Times, stim_Len, @SR_distanceB);
        
        %% Obtain Bird and distance Confusion matrices with MI (Rate distances)
        % add up per line all bins that correspond to the same ID
        ConfMatTempBirdSR=nan(nFiles,NbirdID);
        ConfMatTempDistSR=nan(nFiles,NDist);
        for row=1:nFiles
            for col=1:NbirdID
                colInd = strcmp(verif(:,1),birdID(col));
                ConfMatTempBirdSR(row,col) = sum(ConfMatSR(row,colInd));
            end
            for col=1:NDist
                colInd = strcmp(verif(:,3),StimType(col));
                ConfMatTempDistSR(row,col) = sum(ConfMatSR(row,colInd));
            end
        end
        % Now add up per colon all lines that correspond to the same ID
        % or same distance
        ConfMatBirdSR=nan(NbirdID,NbirdID);
        ConfMatDistSR=nan(NDist,NDist);
        for col=1:NbirdID
            for row=1:NbirdID
                rowInd = strcmp(verif(:,1),birdID(row));
                ConfMatBirdSR(row,col)= sum(ConfMatTempBirdSR(rowInd,col));
            end
        end
        for col=1:NDist
            for row=1:NDist
                rowInd = strcmp(verif(:,3),StimType(row));
                ConfMatDistSR(row,col)= sum(ConfMatTempDistSR(rowInd,col));
            end
        end
        mi_confSR_Bird = info_matrix(ConfMatBirdSR);
        mi_confSR_Dist = info_matrix(ConfMatDistSR);
        
        %% Calculate the ideal oberserver with VR_distanceB for each window size
        for is=1:ns % ns : window size
            
            % Calculate the Van Rossum distances between stimuli and
            % construc a confusion matrix
            [pc, mi_conf, zdT, pzdT, mi_zdT, mi_zdT_nc, zdS, pzdS, mi_zdS, mi_zdS_nc, ConfMat,~,VRDistance_Stims{is},VR_ResponseStims{is}] = info_distanceB_slvg(nFiles, spike_Times, stim_Len, @VR_distanceB, wind(is));
            
            
            percorrectB(is) = pc;
            mi_confusionB(is) = mi_conf;
            zdistTB(is) = zdT;
            pzdistTB(is) = pzdT;
            mizdistTB(is) = mi_zdT;
            mizdistncTB(is) = mi_zdT_nc;
            zdistSB(is) = zdS;
            pzdistSB(is) = pzdS;
            mizdistSB(is) = mi_zdS;
            mizdistncSB(is) = mi_zdS_nc;
            megaMat_bird{is} = ConfMat;
            %added a variable to save the initial confusion matrix (ConfMat) for each window size at 2m here the code only save the last one calculated (so largest window one)
            
            %% Obtain Bird and distance Confusion matrices with MI (VR distances)
            % add up per line all bins that correspond to the same ID
            ConfMatTempBirdVR=nan(nFiles,NbirdID);
            ConfMatTempDistVR=nan(nFiles,NDist);
            for row=1:nFiles
                for col=1:NbirdID
                    colInd = strcmp(verif(:,1),birdID(col));
                    ConfMatTempBirdVR(row,col) = sum(ConfMat(row,colInd));
                end
                for col=1:NDist
                    colInd = strcmp(verif(:,3),StimType(col));
                    ConfMatTempDistVR(row,col) = sum(ConfMat(row,colInd));
                end
            end
            % Now add up per colon all lines that correspond to the same ID
            % or same distance
            ConfMatBirdVR{is}=nan(NbirdID,NbirdID);
            ConfMatDistVR{is}=nan(NDist,NDist);
            for col=1:NbirdID
                for row=1:NbirdID
                    rowInd = strcmp(verif(:,1),birdID(row));
                    ConfMatBirdVR{is}(row,col)= sum(ConfMatTempBirdVR(rowInd,col));
                end
            end
            for col=1:NDist
                for row=1:NDist
                    rowInd = strcmp(verif(:,3),StimType(row));
                    ConfMatDistVR{is}(row,col)= sum(ConfMatTempDistVR(rowInd,col));
                end
            end
            mi_confVR_Bird(is) = info_matrix(ConfMatBirdVR{is});
            mi_confVR_Dist(is) = info_matrix(ConfMatDistVR{is});
            
% Old Code from Solveig
%             %% add up all distances pertaining to one individual in ConfMat (nfiles*nfiles) to obtain a "simple" individual-based ConfMat (NbirdID*NbirdID)
%             ConfMatTemp=zeros(NbirdID);
%             clear BirdMat
%             for row=1:NbirdID
%                 for col=1:NbirdID
%                     BirdMat=ConfMat((row-1)*incr+1:(row-1)*incr+incr, (col-1)*incr+1:(col-1)*incr+incr);
%                     ConfMatTemp(row,col) = sum(sum(BirdMat));
%                 end
%             end
%             if (sum(sum(ConfMatTemp))<0.9999 || sum(sum(ConfMatTemp))>1.0001) % the sum is not always completely equal to 1 for some reason (= to "1.0000"...)
%                 error('ERROR: The sum of probabilities in the confusion matrix is not equal to 1 in ConfMatTemp!');
%             end
%             mi_confTemp = info_matrix(ConfMatTemp);
%             mi_conf_perBird(is) = mi_confTemp;   % mutual information for the Confusion Matrix calculated per bird (per individual to discriminate)
%             ConfMat_bird{is} = ConfMatTemp;      % Confusion Matrix calculated per bird
%             %percorrect_perBird(is) = sum(diag(ConfMat_bird{is}))*100;
            mi_conf_perBird(is) = mi_confVR_Bird(is); % We keep this to make sure subseuqent code using these don't break!
            ConfMat_bird{is} = ConfMatBirdVR{is};
        end
        
        % Finding the winSize that gives best results (i.e., highest mutual information) at 2m, in order to apply it everywhere else
        %         index = find(mi_conf_perBird==max(mi_conf_perBird));
        %         if length(index)>1
        %             fprintf('WARNING! There is more than one max value for the best WinSize for info_distance_B. Taking the first value.\n');
        %             index = index(1);
        %         end
        
        % To keep track of the names of different individuals (right now for Con, 2m only), we create a cell of size (1, ns) or we won't be able to save it in the output
        % the variable to keep is stim_IDR{1}
        stimID = cell(1,ns);
        for i = 1:ns
            stimID{i} = Stim_IDR{1};
        end
        megaMat_type = 'per individual';
        Verif = {verif};
        
        % Output structure. Important values for later analysis are: stimrate, avgzscore, mean_dprime, baseline_dprime, mi_conf_perBird.
        output = struct();
        output.StimType = ConSyn(qq);
        output.megaMat_type = megaMat_type;
        output.nfiles = nFiles;
        output.winsize=wind;
        output.mi_conf_perBird = mi_conf_perBird;
        output.Percorrect_perBird = percorrect_perBird;
        output.ConfMatBird = ConfMat_bird;
        output.megaMat_bird = megaMat_bird;
        output.StimID_2m = stimID;
        output.verif = Verif;
        output.StimDurationms = stim_Len;
        output.VRDistanceStims = VRDistance_Stims;
        output.SRDistanceStims = SRDistance_Stims;
        output.SR_ResponseStims = SR_ResponseStims;
        output.VR_ResponseStims = VR_ResponseStims;
        output.spike_Rate = spike_Rate;
        output.mi_confusionB = mi_confusionB;
        output.mi_confusionB_SR = mi_confSR;
        output.ConfMatBirdSR = ConfMatBirdSR;
        output.mi_confusionB_SR_Bird = mi_confSR_Bird;
        output.ConfMatDistSR = ConfMatDistSR;
        output.mi_confusionB_SR_Dist = mi_confSR_Dist;
        output.ConfMatBirdVR = ConfMatBirdVR;
        output.mi_confusionB_VR_Bird = mi_confVR_Bird;
        output.ConfMatDistVR = ConfMatDistVR;
        output.mi_confusionB_VR_Dist = mi_confVR_Dist;
        
        
        
        %megaData{1} = output; % data of megaMatrix grouped per individual
        dataMega{qq} = output;
        
        %% for 'Con' calls, we calculate the mutual info for all window Sizes, regrouping stims per distance for all birds
        
        %         clear spike_Times stim_Len verif mi_conf_perBird
        %         percorrectB = zeros(1,ns);
        %         mi_confusionB = zeros(1,ns);
        %         zdistTB = zeros(1,ns);              % Uses std for other for both self and other in KL Divergence calc
        %         pzdistTB = zeros(1,ns);
        %         mizdistTB = zeros(1,ns);
        %         mizdistncTB = zeros(1,ns);
        %         zdistSB = zeros(1,ns);              % Recenters "other" distrubtions to zero to minimize additive variance across songs
        %         pzdistSB = zeros(1,ns);             % and collapses distances across all trials of a single stimulus
        %         mizdistSB = zeros(1,ns);
        %         mizdistncSB = zeros(1,ns);
        %         megaMat_dist = cell(1,ns);
        %         ConfMat_dist = cell(1,ns);
        %         mi_conf_perDist = zeros(1,ns);
        %         percorrect_perDist = zeros(1,ns);
        %
        %         % organise the megaMatrix per distance
        %         incr = ncalls*NbirdID; % nb of regrouped calls per distance
        %         spike_Times = cell(1, NselectorInd*NDist);
        %         stim_Len = zeros(1, NselectorInd*NDist);
        %         verif = cell(NselectorInd*NDist, 3);
        %         for ii = 1:NDist
        %             Trains_temp = cell(1, incr);
        %             Len_temp = zeros(1, incr);
        %             verif_temp = cell(incr, 3);
        %             for jj = 1:NbirdID
        %                 Trains_temp((jj-1)*ncalls+1:(jj-1)*ncalls+ncalls) = SpikeTrainsR{ii}((jj-1)*ncalls+1:(jj-1)*ncalls+ncalls);
        %                 Len_temp((jj-1)*ncalls+1:(jj-1)*ncalls+ncalls) = Sections_lenR{ii}((jj-1)*ncalls+1:(jj-1)*ncalls+ncalls);
        %                 verif_temp((jj-1)*ncalls+1:(jj-1)*ncalls+ncalls, :) = Stim_IDR{ii}((jj-1)*ncalls+1:(jj-1)*ncalls+ncalls, :);
        %             end
        %             spike_Times((ii-1)*incr+1:(ii-1)*incr+(incr)) = Trains_temp;
        %             stim_Len((ii-1)*incr+1 : (ii-1)*incr+(incr)) = Len_temp;
        %             verif((ii-1)*incr+1:(ii-1)*incr+(incr), :) = verif_temp;  % pour vérification de la classif
        %         end
        %         nFiles=length(spike_Times);
        %
        %         % Calculate the ideal oberserver with VR_distanceB
        %         for is=1:ns % ns : window size
        %
        %             % Reinsert this if using the confusion matrix in version B
        %             [pc, mi_conf, zdT, pzdT, mi_zdT, mi_zdT_nc, zdS, pzdS, mi_zdS, mi_zdS_nc, ConfMat] = info_distanceB_slvg(nFiles, spike_Times, stim_Len, @VR_distanceB, wind(is));
        %
        %             percorrectB(is) = pc;
        %             mi_confusionB(is) = mi_conf;
        %             zdistTB(is) = zdT;
        %             pzdistTB(is) = pzdT;
        %             mizdistTB(is) = mi_zdT;
        %             mizdistncTB(is) = mi_zdT_nc;
        %             zdistSB(is) = zdS;
        %             pzdistSB(is) = pzdS;
        %             mizdistSB(is) = mi_zdS;
        %             mizdistncSB(is) = mi_zdS_nc;
        %             megaMat_dist{is} = ConfMat; %added a variable to save the initial confusion matrix (ConfMat) for each window size at 2m here the code only save the last one calculated (so largest window one)
        %
        %             % add up all distances pertaining to one individual in ConfMat (nfiles*nfiles) to obtain a "simple" individual-based ConfMat (NbirdID*NbirdID)
        %             ConfMatTemp=zeros(NDist);
        %             clear BirdMat
        %             for row=1:NDist
        %                 for col=1:NDist
        %                     BirdMat=ConfMat((row-1)*incr+1:(row-1)*incr+incr, (col-1)*incr+1:(col-1)*incr+incr);
        %                     ConfMatTemp(row,col) = sum(sum(BirdMat));
        %                 end
        %             end
        %             if (sum(sum(ConfMatTemp))<0.9999 || sum(sum(ConfMatTemp))>1.0001) % the sum is not always completely equal to 1 for some reason (= to "1.0000"...)
        %                 error('ERROR: The sum of probabilities in the confusion matrix is not equal to 1 in ConfMatTemp!');
        %             end
        %             mi_confTemp = info_matrix(ConfMatTemp);
        %             mi_conf_perDist(is) = mi_confTemp;   % mutual information for the Confusion Matrix calculated per bird (per individual to discriminate)
        %             ConfMat_dist{is} = ConfMatTemp;      % Confusion Matrix calculated per bird
        %             percorrect_perDist(is) = sum(diag(ConfMat_bird{is}))*100;
        %         end
        %
        %         % Finding the winSize that gives best results (i.e., highest mutual information) at 2m, in order to apply it everywhere else
        % %         index = find(mi_conf_perBird==max(mi_conf_perBird));
        % %         if length(index)>1
        % %             fprintf('WARNING! There is more than one max value for the best WinSize for info_distance_B. Taking the first value.\n');
        % %             index = index(1);
        % %         end
        %
        %         % To keep track of the names of different individuals (right now for Con, 2m only), we create a cell of size (1, ns) or we won't be able to save it in the output
        %         % the variable to keep is stim_IDR{1}
        %         stimID = cell(1,ns);
        %         for i = 1:ns
        %             stimID{i} = Stim_IDR{1};
        %         end
        %         megaMat_type = 'per distance';
        %         Verif = {verif};
        %
        %         % Output structure. Important values for later analysis are: stimrate, avgzscore, mean_dprime, baseline_dprime, mi_conf_perBird.
        %         output = struct('StimType', ConSyn(qq), 'megaMat_type', megaMat_type, 'nfiles', nFiles, 'winsize', wind, ...
        %             'mi_conf_perDist', mi_conf_perDist, 'Percorrect_perDist', percorrect_perDist, 'ConfMatDist', ConfMat_dist, ...
        %             'megaMat_dist', megaMat_dist, 'StimID_2m', stimID, 'verif', Verif);
        %
        %         megaData{2} = output; % data of megaMatrix grouped per individual
        
        %dataMega{qq} = output;
    end
end
%%
[pathstr,h5Name,ext] = fileparts(h5Path);
fileName = strcat('MegaMat_', Birdname, '_', h5Name, '.mat');
save(fullfile(outputFolder, fileName), 'dataMega');

return






