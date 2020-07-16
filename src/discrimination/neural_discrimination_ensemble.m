function outputs = neural_discrimination_ensemble(neuron_id, stimtype, maxrep, stim_names_test) 
% Calculates 3 measures of neural discrimination, the within dprime, the
% percent correct of an ideal oberver based on VR method and the
% information based on Gamma model of spiking neuron
% read all the spike arrival times in a folder
% Input: n_neurons = number of neurons
%        neuron_id = a structure that has the length of n_neurons and
%        elements birdname, brainregion, cellname

n_neurons = length(neuron_id);              % number of neurons
n_files = zeros(1, n_neurons);
spike_times = cell(1, n_neurons);
stim_len = cell(1, n_neurons);
stim_name = cell(1, n_neurons);

% Read all the data
for i_neuron=1:n_neurons
% Read input data from data base
    [n_files(i_neuron) spike_times{i_neuron} stim_len{i_neuron} stim_name{i_neuron}] = ...
                    read_all_spikes(neuron_id(i_neuron).birdname, neuron_id(i_neuron).brainregion, neuron_id(i_neuron).cellname, stimtype);
end

% Find common stimuli
n_good_stims = 0;

% Check if stimuli are requested and if not choose the stims from first
% neuron 
if isempty(stim_names_test)
   i_neuron = 1;
   for i_file=1:n_files(i_neuron)
    stim_names_test{i_file} = stim_name{i_neuron}(i_file);
   end
end
n_test = length(stim_names_test);
   
for i_file=1:n_test;
    stim_hits = 0;
    stim_ind_test = zeros(1, n_neurons);
    for i2_neuron=1:n_neurons
        for i2_file=1:n_files(i2_neuron)
            if strcmp(stim_names_test{i_file}, stim_name{i2_neuron}(i2_file))
                stim_hits = stim_hits+1;
                stim_ind_test(i2_neuron) = i2_file;
                break;
            end
        end
    end
    if stim_hits == n_neurons
        n_good_stims = n_good_stims + 1;
        stim_ind(n_good_stims, :) = stim_ind_test;
        good_stim_name{n_good_stims} = stim_names_test{i_file};
    end
end

if n_good_stims == 0 
    fprintf(1, 'Error in neural_discrimination_ensemble. No common stimulus:\n');
    for i2_neuron=1:n_neurons
        fprintf(1, 'Neuron %d: %s %s has %d stims: \n', i2_neuron, neuron_id(i2_neuron).birdname, neuron_id(i2_neuron).cellname, n_files(i2_neuron));
        for i2_file=1:n_files(i2_neuron)
            fprintf(1, '%s ', stim_name{i2_neuron}{i2_file});
        end
        fprintf(1, '\n');
    end
    outputs = [];
    return;
end

stim_len_ensemble = stim_len{1}(stim_ind(1:n_good_stims, 1));

        
fprintf(1,'n_neurons = %d n_good_stims = %d\n', n_neurons, n_good_stims);
for i_neuron = 1:n_neurons
    fprintf( 1, '\tBird = %s Cell Name = %s Stims ID:\n\t', neuron_id(i_neuron).birdname, neuron_id(i_neuron).cellname);
    for i_good_stim=1:n_good_stims
        fprintf( 1,'%d ', stim_ind(i_good_stim, i_neuron));
        nrep = length(spike_times{i_neuron}{stim_ind(i_good_stim, i_neuron)});
        if (nrep < maxrep)
            maxrep = nrep;
        end
    end
    fprintf(1, '\n');
end
fprintf(1, 'Maximun number of trials common to all cells %d\n', maxrep); 
 
%  Allocate space for good stims in the correct order
spike_times_ensemble = cell(1, n_good_stims);
spike_times_ensemble_for_rate = cell(1, n_good_stims);
for i_good_stim=1:n_good_stims
    spike_times_ensemble{i_good_stim} = cell(1, maxrep);
    for irep=1:maxrep
        spike_times_ensemble{i_good_stim}{irep} = cell(1, n_neurons);
        for ineuron=1:n_neurons
            spike_times_ensemble{i_good_stim}{irep}{ineuron} = spike_times{ineuron}{stim_ind(i_good_stim, ineuron)}{irep};
            spike_times_ensemble_for_rate{ineuron}{i_good_stim}{irep} = spike_times{ineuron}{stim_ind(i_good_stim, ineuron)}{irep};            
        end
    end
end


% Window sizes for calculations that depend on window length
winSize = [5];    % Window size for static Info calculation and for ideal observer
% winSize = [50];    % Window size for static Info calculation and for ideal observer
ns = length(winSize);


%***Ideal Observer Calculations (by Trial Template)***
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
confusions = cell(1,ns);



if n_good_stims
    
    for i_neuron = 1:n_neurons
       [stimrate(i_neuron) stdstimrate(i_neuron) backrate(i_neuron) stdbackrate(i_neuron) avgzscore(i_neuron) stdzscore(i_neuron) avgdprime(i_neuron) stddprime(i_neuron)] = ...
                dprime_within(n_good_stims, spike_times_ensemble_for_rate{i_neuron}, stim_len_ensemble);
    end
    
    % With info_distance B calling VR_multidistance
    for is=1:ns     
        [pc, mi_conf, zdT, pzdT, mi_zdT, mi_zdT_nc, zdS, pzdS, mi_zdS, mi_zdS_nc, confusion_matrix] = info_distanceB(n_good_stims, spike_times_ensemble, stim_len_ensemble, @VR_multidistance, winSize(is));
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
        confusions{is} = confusion_matrix;
    end
    


end


outputs = struct('nfiles', n_good_stims, 'winsize', winSize, 'maxrep', maxrep, 'n_good_stims', n_good_stims, 'n_neurons', n_neurons, 'stimrate', stimrate, ...
    'percorrectB', percorrectB, 'mi_confusionB', mi_confusionB, 'zdistTB', zdistTB, 'pzdistTB', pzdistTB, 'mizdistTB', mizdistTB, 'mizdistncTB', mizdistncTB, ...
                                                                'zdistSB', zdistSB, 'pzdistSB', pzdistSB, 'mizdistSB', mizdistSB, 'mizdistncSB', mizdistncSB, 'confmat', confusions);
outputs.stim_names = good_stim_name;                                                            

return





