% Calculate MI of ensembles and single cells
clear

n_ensemble = 2:10;
n_trials = 100;

for i = 1:length(n_ensemble)
    for trial = 1:n_trials
%                 [n, output_ensemble, output_single]=...
%                      Normal_L_discrimination_ensemble(n_ensemble(i));
        
        [n, output_ensemble, output_single]=...
            wnL_discrimination_ensemble(n_ensemble(i));
        
        % Create a structure containing: n_ensemble, ensemble MI, and
        % single MI
        info.n_ensemble = n;
        if isempty(output_ensemble)
            ensemble_info{i}{trial} = [];
        else
            info.ensembleMI = output_ensemble.mizdistSB;
            info.singleMI = [output_single.mizdistSB];
            info.rate = output_ensemble.stimrate;
            
            % Store structure in cell array: ensemble_info
            ensemble_info{i}{trial} =  info;
        end
    end
end

save ensemble_info_t100Rate ensemble_info
