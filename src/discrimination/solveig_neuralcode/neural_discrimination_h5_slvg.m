function data = neural_discrimination_h5_slvgFinal(h5Path, Birdname)
%[nfiles, stimrate, stdstimrate, backrate, stdbackrate, avgzscore, stdzscore, avgdprime, stddprime, percorrect, avgrank, info, noiseinfo, rateinfo, info1, info2, info2Interval]
% Calculates 3 measures of neural discrimination, the within dprime, the
% percent correct of an ideal oberver based on VR method and the
% information based on Gamma model of spiking neuron
% read all the spike arrival times in a folder
% stimType : 'Mask1' etc. Was in the function input in previous versions of this script, not anymore in the 'Solveig' version (see l.23).

%Birdname = 'GraBla1602F'; % example...
%h5Path = 'Site1_R2601_e10_s0_ss1.h5';

%% Read input data from data base
unit = read_unit_h5file(h5Path, 'r'); % 'r' to read only, 'a' to add fields

%% Select the protocol (SelX, Callx, Maskx, STRFx...)
%number of different protocols run for that unit, identify if the one
%you are asking for is existing in this unit and selecting this one in
%responses
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
nresp = size(responses,2);

%% Construct the input/output we need for all the following code
% This is the number of sound files played
nfiles = length(responses);

%This is the vector of the length of the different stims
stim_len = zeros (nfiles,1);
stim_class = cell (nfiles,1);
birdid = cell (nfiles,1);
distance = cell (nfiles,1);

for nf=1:nfiles
    stim_len(nf) = 1000*responses{nf}.stim_duration;%here stim_len is in ms
    stim_class{nf} = responses{nf}.stim_class;
    birdid{nf} = responses{nf}.birdid;
    distance{nf} = responses{nf}.distance;
end

% This is the cell array of the time of spike arrival
spike_times = cell(1, nfiles);
for nf = 1:nfiles
    Trials = responses{nf}.trials;
    nt = length(Trials);
    spike_times{nf} = cell(1, nt);
    for it = 1:nt
        trial = Trials{it};
        spike_times{nf}{it} = trial.spikeTimes.*1000.0;   % spike_times are in ms in neural_discrimination and in s in hdf5
    end
end

% This is various stuff we'll need for the final outputs
ConSyn = {'Con', 'Syn'};
data = cell(1, 2);      % Final data cell (Con and Syn)
    
%% Set up parameters for all the calculations    
% Window sizes for calculations that depend on window length
winSize = [1 3 5 10 30 50 100];    % Window size for static Info calculation and for ideal observer
% We will run the code for 2m for 'Con' calls with all winSizes, choose the best window and apply it for all distances and 'Syn' calls.
ns = length(winSize);

% Reorder the distance vector from 2 m to 256 m
dist=unique(distance); % !!!!! wrong order : 128, 16, 2, 256, 64 m.
for i=1:length(dist)
di(i) = str2num(dist{i});           % change string to number...
end
di = sort(di);                          % reorder the distances...
for i=1:length(dist)
StimType{i} = num2str(di(i));     % put it back into a cell of strings
end
NstimType=length(StimType); % = number of distances


%% Calculations!!

if nfiles
    for qq = [1 2] % 1: Conspecific calls, 2 : Synthetic calls
        
        if qq==1 % Compare bird identities depending on distance, for Conspecific calls
            NDist = NstimType; % all distances for 'Con' calls
            dataPerDist = cell(1, NDist);
            PropDist = {'2m', '16m', '64m', '128m', '256m'};
            
            % Create the data set for each distance separately, using only Conspecific calls
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
                Stim_ID_local = cell(NselectorInd, 3);
                for NSI = 1:NselectorInd 
                    SI=selectorInd(NSI);
                    SpikeTrainslocal{NSI} = spike_times{SI};
                    Sections_lenlocal(NSI)= stim_len(SI);
                    Stim_ID_local{NSI, 1} = responses{SI}.birdid;
                    Stim_ID_local{NSI, 2} = responses{SI}.callnum;
                    Stim_ID_local{NSI, 3} = responses{SI}.distance;
                end
                SpikeTrains{vt}=SpikeTrainslocal; % changer le nom ou faire une méta-structure si on fait tourner la boucle Con-Syn...
                Sections_len{vt}=Sections_lenlocal;
                Stim_ID{vt}=Stim_ID_local;
            end

            % Reorganise the cells to have all calls from the same birdid in order
            SpikeTrainsR=cell(NDist,1);
            Sections_lenR=cell(NDist,1);
            Stim_IDR=cell(NDist, 1);
            for vt=1:NDist
                birdID=unique(Stim_ID{vt}(:,1));
                NbirdID=length(birdID);
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
                        Sections_lenR{vt}(incr+ind) = Sections_len{vt}(selectorId(ind));
                    end
                end
            end
    
    
        elseif qq==2 % Compare bird identities depending on distance, for Synthetic calls
            NDist = NstimType-1; % All distances except 2m for 'Syn' calls
            dataPerDist = cell(1, NDist);
            PropDist = {'16m', '64m', '128m', '256m'};
            
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

            % Reorganise the cells to have al calls from the same birdid in order
            SpikeTrainsR=cell(NDist,1);
            Sections_lenR=cell(NDist,1);
            Stim_IDR=cell(NDist, 1);
            for vt=1:NDist
                birdID=unique(Stim_ID{vt}(:,1));
                NbirdID=length(birdID);
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
                        Sections_lenR{vt}(incr+ind) = Sections_len{vt}(selectorId(ind));
                    end
                end
            end
        end
        
                %% Initialize all return values to zero
                %***DPrime Calculations***
                stimrate = 0;
                stdstimrate = 0;
                backrate = 0;
                stdbackrate = 0;
                avgzscore = 0;
                stdzscore = 0;
                avgdprime = 0;
                stddprime = 0;
                dprime_mat = [];
                
        %% ********Start of the calculations!*********** 
                
        for Dist=1:NDist
            
            if (qq==1 & Dist==1) 
            %% for 2m and 'Con' calls, we calculate the mutual info for all window Sizes
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
                ConfMatBird = cell(1,ns);
                mi_conf_perBird = zeros(1,ns);
                ConfMatCall = cell(1,ns);
                
                spike_Times=SpikeTrainsR{Dist};
                stim_Len=Sections_lenR{Dist};
                nFiles=length(spike_Times);
                % Calculate dprime within
                [stimrate stdstimrate backrate stdbackrate avgzscore stdzscore avgdprime stddprime dprime_mat] = dprime_within_slvg(nFiles, spike_Times, stim_Len);
                % !! "avgdprime" is not the correct value here, as it is calculated without taking into account the bird individuals.

                % Calculate the mean(dprime) for each individual in the dprime_matrix, then the overall average dprime (that will be called mean_dprime) for all 
                % 2 by 2 comparisons (not taking the diagonal into account)
                d_mat_temp = zeros(NbirdID);
                clear BirdMat
                ncalls=nFiles/NbirdID;
                if (ncalls ~= round(ncalls)) % check that ncalls is an integer number
                    error('ERROR: The number of calls per individual is not a multiple of the total number of calls!');
                end
                for row=1:NbirdID
                    for col=1:NbirdID
                        BirdMat = dprime_mat((row-1)*ncalls+1:(row-1)*ncalls+ncalls, (col-1)*ncalls+1:(col-1)*ncalls+ncalls);
                        if row==col
                            d_mat_temp(row, col) = (sum(sum(BirdMat)))/(ncalls*(ncalls-1)); % mean(dprime) per bird, leaving out the diagonal of zeros
                        else
                            d_mat_temp(row, col) = (sum(sum(BirdMat)))/(ncalls*ncalls); % mean(dprime) for different birds
                        end
                    end
                end
                mean_dprime = (sum(sum(d_mat_temp))-sum(diag(d_mat_temp)))/(NbirdID*(NbirdID-1));
                % calculating the std deviation for this value... (sheesh)
                dpri = d_mat_temp;
                for i = 1:NbirdID
                    dpri(i,i) = -1; % give negative value to numbers in diagonal
                end
                d_mat = reshape(dpri, NbirdID*NbirdID, 1);
                d_mat(d_mat<0) = []; % remove negative values
                std_dprime = std(d_mat);

                % calculating the mean dprime for same bird comparison (mean of the values in the diagonal of the matrix)
                baseline_dprime = sum(diag(d_mat_temp))/NbirdID;
                std_base_dprime = std(diag(d_mat_temp));

                
                % Calculate the ideal oberserver with VR_distanceB
                for is=1:ns % ns : window size

                    % Reinsert this if using the confusion matrix in version B
                    [pc, mi_conf, zdT, pzdT, mi_zdT, mi_zdT_nc, zdS, pzdS, mi_zdS, mi_zdS_nc, ConfMat] = info_distanceB_slvg(nFiles, spike_Times, stim_Len, @VR_distanceB, winSize(is));

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
                    ConfMatCall{is} = ConfMat; %added a variable to save the initial confusion matrix (ConfMat) for each window size at 2m here the code only save the last one calculated (so largest window one)

                    % add up all distances pertaining to one individual in ConfMat (nfiles*nfiles) to obtain a "simple" individual-based ConfMat (NbirdID*NbirdID)
                    ConfMatTemp=zeros(NbirdID);
                    clear BirdMat
                    for row=1:NbirdID
                        for col=1:NbirdID
                            BirdMat=ConfMat((row-1)*ncalls+1:(row-1)*ncalls+ncalls, (col-1)*ncalls+1:(col-1)*ncalls+ncalls);
                            ConfMatTemp(row,col) = sum(sum(BirdMat));
                        end
                    end
                    if (sum(sum(ConfMatTemp))<0.9999 || sum(sum(ConfMatTemp))>1.0001) % the sum is not always completely equal to 1 for some reason (= to "1.0000"...)
                        error('ERROR: The sum of probabilities in the confusion matrix is not equal to 1 in ConfMatTemp!');
                    end
                    mi_confTemp = info_matrix(ConfMatTemp);
                    mi_conf_perBird(is) = mi_confTemp;  % mutual information for the Confusion Matrix calculated per bird (per individual to discriminate)
                    ConfMatBird{is} = ConfMatTemp;      % Confusion Matrix calculated per bird
                end
                
                % Finding the winSize that gives best results (i.e., highest mutual information) at 2m, in order to apply it everywhere else
                index = find(mi_conf_perBird==max(mi_conf_perBird));
                if length(index)>1
                    fprintf('WARNING! There is more than one max value for the best WinSize for info_distance_B. Taking the first value.\n');
                    index = index(1);
                end
                
                % To keep track of the names of different individuals (right now for Con, 2m only), we create a cell of size (1, ns) or we won't be able to save it in the output
                % the variable to keep is stim_IDR{1}
                stimID = cell(1,ns);
                for i = 1:ns
                    stimID{i} = Stim_IDR{1};
                end
                
                % Output structure. Important values for later analysis are: stimrate, avgzscore, mean_dprime, baseline_dprime, mi_conf_perBird.
                output = struct('StimType', ConSyn(qq), 'Distance', PropDist(Dist), 'nfiles', nFiles, 'winsize', winSize, 'ChosenWinSize', winSize(index), 'stimrate', stimrate, 'stdstimrate', stdstimrate, ...
                    'backrate', backrate, 'stdbackrate', stdbackrate, 'avgzscore', avgzscore, 'stdzscore', stdzscore, 'mean_dprime', mean_dprime, 'std_dprime', std_dprime, ...
                    'baseline_dprime', baseline_dprime, 'std_base_dprime', std_base_dprime, 'mi_conf_perBird', mi_conf_perBird, 'ConfMatBird', ConfMatBird, 'percorrect', percorrectB, ...
                    'mi_confusion', mi_confusionB, 'zdistTB', zdistTB, 'pzdistTB', pzdistTB, 'mizdistTB', mizdistTB, 'mizdistncTB', mizdistncTB, 'zdistSB', zdistSB, 'pzdistSB', pzdistSB, ...
                    'mizdistSB', mizdistSB, 'mizdistncSB', mizdistncSB, 'ConfMatTot', ConfMatCall, 'StimID_2m', stimID);
                
                dataPerDist{Dist} = output;
                
            else
                %% for the rest of distances and for Syn calls, we use the winSize that gave the best results at 2m.
                
                spike_Times=SpikeTrainsR{Dist};
                stim_Len=Sections_lenR{Dist};
                nFiles=length(spike_Times);
                % Calculate dprime within
                [stimrate stdstimrate backrate stdbackrate avgzscore stdzscore avgdprime stddprime dprime_mat] = dprime_within_slvg(nFiles, spike_Times, stim_Len);
                % !! "avgdprime" is not the correct value here, as it is calculated without taking into account the bird individuals.
                
                % Calculate the mean(dprime) for each individual in the dprime_matrix, then the overall average dprime (that will be called mean_dprime) for all
                % 2 by 2 comparisons (not taking the diagonal into account)
                d_mat_temp = zeros(NbirdID);
                clear BirdMat
                ncalls=nFiles/NbirdID;
                if (ncalls ~= round(ncalls)) % check that ncalls is an integer number
                    error('ERROR: The number of calls per individual is not a multiple of the total number of calls!');
                end
                for row=1:NbirdID
                    for col=1:NbirdID
                        BirdMat = dprime_mat((row-1)*ncalls+1:(row-1)*ncalls+ncalls, (col-1)*ncalls+1:(col-1)*ncalls+ncalls);
                        if row==col
                            d_mat_temp(row, col) = (sum(sum(BirdMat)))/(ncalls*(ncalls-1)); % mean(dprime) per bird, leaving out the diagonal of zeros
                        else
                            d_mat_temp(row, col) = (sum(sum(BirdMat)))/(ncalls*ncalls); % mean(dprime) for different birds
                        end
                    end
                end
                mean_dprime = (sum(sum(d_mat_temp))-sum(diag(d_mat_temp)))/(NbirdID*(NbirdID-1));
                % calculating the std deviation for this value... (sheesh)
                dpri = d_mat_temp;
                for i = 1:NbirdID
                    dpri(i,i) = -1; % give negative value to numbers in diagonal
                end
                d_mat = reshape(dpri, NbirdID*NbirdID, 1);
                d_mat(d_mat<0) = []; % remove negative values
                std_dprime = std(d_mat);
                
                % calculating the mean dprime for same bird comparison (mean of the values in the diagonal of the matrix)
                baseline_dprime = sum(diag(d_mat_temp))/NbirdID;
                std_base_dprime = std(diag(d_mat_temp));
                
                
                % Calculate the ideal oberserver with VR_distanceB
                for is = index % window size, chosen at 2m
                    
                    % Reinsert this if using the confusion matrix in version B
                    [pc, mi_conf, zdT, pzdT, mi_zdT, mi_zdT_nc, zdS, pzdS, mi_zdS, mi_zdS_nc, ConfMatCall] = info_distanceB_slvg(nFiles, spike_Times, stim_Len, @VR_distanceB, winSize(is));
                    
                    percorrectB(is) = pc;
                    mi_confusionB(is) = mi_conf;
                    
                    zdistTB = zdT;
                    pzdistTB = pzdT;
                    mizdistTB = mi_zdT;
                    mizdistncTB = mi_zdT_nc;
                    zdistSB = zdS;
                    pzdistSB = pzdS;
                    mizdistSB = mi_zdS;
                    mizdistncSB = mi_zdS_nc;

                    % add up all distances pertaining to one individual in ConfMat (nfiles*nfiles) to obtain a "simple" individual-based ConfMat (NbirdID*NbirdID)
                    ConfMatTemp=zeros(NbirdID);
                    clear BirdMat
                    for row=1:NbirdID
                        for col=1:NbirdID
                            BirdMat=ConfMatCall((row-1)*ncalls+1:(row-1)*ncalls+ncalls, (col-1)*ncalls+1:(col-1)*ncalls+ncalls);
                            ConfMatTemp(row,col) = sum(sum(BirdMat));
                        end
                    end
                    if (sum(sum(ConfMatTemp))<0.9999 || sum(sum(ConfMatTemp))>1.0001) % the sum is not always completely equal to 1 for some reason (= to "1.0000"...)
                        error('ERROR: The sum of probabilities in the confusion matrix is not equal to 1 in ConfMatTemp!');
                    end
                    mi_confTemp = info_matrix(ConfMatTemp);
                    mi_conf_perBird = mi_confTemp;
                    ConfMatBird = ConfMatTemp;
                end
                
            % Output structure. Important values for later analysis are: stimrate, avgzscore, mean_dprime, baseline_dprime, mi_conf_perBird.
            output = struct('StimType', ConSyn(qq), 'Distance', PropDist(Dist), 'nfiles', nFiles, 'winsize', winSize, 'stimrate', stimrate, 'stdstimrate', stdstimrate, 'backrate', backrate, 'stdbackrate', stdbackrate, ...
            'avgzscore', avgzscore, 'stdzscore', stdzscore, 'mean_dprime', mean_dprime, 'std_dprime', std_dprime, 'baseline_dprime', baseline_dprime, 'std_base_dprime', std_base_dprime, ...
            'mi_conf_perBird', mi_conf_perBird, 'ConfMatBird', ConfMatBird, 'percorrect', percorrectB, 'mi_confusion', mi_confusionB, 'zdistTB', zdistTB, 'pzdistTB', pzdistTB, ...
            'mizdistTB', mizdistTB, 'mizdistncTB', mizdistncTB, 'zdistSB', zdistSB, 'pzdistSB', pzdistSB, 'mizdistSB', mizdistSB, 'mizdistncSB', mizdistncSB, ...
            'ConfMatTot', ConfMatCall);
        
            dataPerDist{Dist} = output;
        
            end
            data{qq} = dataPerDist;
        end
        
    end
end

h5Name = h5Path(1:end-3);
h5Name = strcat('-', h5Name);
fileName = strcat(Birdname, h5Name);
save(fileName, 'data');

return






