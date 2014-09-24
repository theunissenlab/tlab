function [meanSnip, sdSnip] = meanSnip_h5(unit)
% Returns the mean snipet for an h5 file
% unit is the structure returned by read_unit_h5file

ntsnip = 18;       % Lenth of our spike snippets
% Reads the h5 file and calculates the mean snippet

source_directory = unit.source_directory;


% Loop through data structure and get some stats
nfilesTot = 0;
ntrialsTot = 0;
nsnipsTot =0;

nclasses = length(unit.classes);
nfiles = zeros(nclasses,1);
for iclass=1:nclasses
    nfiles(iclass) = length(unit.class_responses.(unit.classes{iclass})); 
    nfilesTot = nfilesTot + nfiles(iclass);     
    for nfi=1:nfiles(iclass)   
        response = unit.class_responses.(unit.classes{iclass}){nfi};    
        ntrials = length(response.trials); 
        ntrialsTot = ntrialsTot + ntrials;
        
        for it=1:ntrials
            trial = response.trials{it};
            spike_id_trials = trial.spikeIds;
            ns = length(spike_id_trials);
            nsnipsTot = nsnipsTot + ns;
        end
    end
end
fprintf(1, 'File statistics:\n\tNumber of Stims %d\n\tNumber of trials %d\n\tNumber of snips %d\n', nfilesTot, ntrialsTot, nsnipsTot);

% Read all the snippets

snips = zeros(nsnipsTot,ntsnip);

cnt=0;
for iclass=1:nclasses
    a = dir([source_directory '/' unit.site '*' unit.classes{iclass} '*waves.f32']);
    fprintf(1, 'Reading and storing snipets from wave file: %s\n',  [source_directory '/' a.name]);
    fsnip = fopen([source_directory '/' a.name], 'r');
    for nfi=1:nfiles(iclass)   
        response = unit.class_responses.(unit.classes{iclass}){nfi};    
        ntrials = length(response.trials); 
        
        for it=1:ntrials
            trial = response.trials{it};
            spike_id_trials = trial.spikeIds;
            ns = length(spike_id_trials);
            for is=1:ns
                cnt = cnt+1;
                fseek(fsnip, (spike_id_trials(is)-1)*4*ntsnip, 'bof');
                snips(cnt,:) = fread(fsnip, ntsnip, 'float32');
                % Note that we must use the full name if we want to store
                % it.
                unit.class_responses.(unit.classes{iclass}){nfi}.trials{it}.spikeWav(is,:) = snips(cnt,:);
                if mod(cnt, 5000)==0
                    fprintf('.')
                end
            end
        end
    end
    fclose(fsnip);
end

fprintf(1, '\n...Done reading snippets\n');

% Calculating the mean snippet
meanSnip = mean(snips,1);
sdSnip = std(snips, 0, 1);

return
