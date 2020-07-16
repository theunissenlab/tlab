
%% Information vs. Ensemble Size
load ensemble_info_t100Rate
n_trials = 100;
n_ensemble = 2:10;
ensembleMI_wn = zeros(length(n_ensemble), n_trials);
redundancy_wn = zeros(length(n_ensemble), n_trials);
sum_single_wn = zeros(length(n_ensemble), n_trials);
ensembleMIspike_wn = zeros(length(n_ensemble), n_trials);
sumspike_single_wn = zeros(length(n_ensemble), n_trials);


for i = 1:length(n_ensemble)
    for trial = 1:n_trials
        %Extract data from cell array
        ensembleMI_wn(i,trial) = ensemble_info{i}{trial}.ensembleMI;
        sum_single_wn(i,trial) = sum(ensemble_info{i}{trial}.singleMI);
        ensembleMIspike_wn(i,trial) = ensembleMI_wn(i,trial)./sum(ensemble_info{i}{trial}.rate);
        sumspike_single_wn(i,trial) = sum(ensemble_info{i}{trial}.singleMI./ensemble_info{i}{trial}.rate);       
        redundancy_wn(i,trial) = (sum_single_wn(i,trial) - ensembleMI_wn(i,trial))/sum_single_wn(i,trial);
    end
end

% load ensemble_info_L
load ensemble_info_L_113 
ensembleMI_control = zeros(length(n_ensemble), n_trials);
redundancy_control = zeros(length(n_ensemble), n_trials);
sum_single_control = zeros(length(n_ensemble), n_trials);
ensembleMIspike_control = zeros(length(n_ensemble), n_trials);
sumspike_single_control = zeros(length(n_ensemble), n_trials);

for i = 1:length(n_ensemble)
    for trial = 1:n_trials
        if isempty(ensemble_info{i}{trial})
            ensembleMI_control(i,trial) = NaN;
            sum_single_control(i,trial) = NaN;
            ensembleMIspike_control(i,trial) = NaN;
            sumspike_single_control(i,trial) = NaN;
            redundancy_control(i,trial) = NaN;
        else          
            %Extract data from cell array
            ensembleMI_control(i,trial) = ensemble_info{i}{trial}.ensembleMI;
            sum_single_control(i,trial) = sum(ensemble_info{i}{trial}.singleMI);
            ensembleMIspike_control(i,trial) = ensembleMI_control(i,trial)./sum(ensemble_info{i}{trial}.rate);
            sumspike_single_control(i,trial) = sum(ensemble_info{i}{trial}.singleMI./ensemble_info{i}{trial}.rate); 
            redundancy_control(i,trial) = (sum_single_control(i,trial) - ensembleMI_control(i,trial))/sum_single_control(i,trial);
        end
    end
end

% Ensemble information plot
average_wn = nanmean(ensembleMI_wn,2);               %Find average information across trials
average_control = nanmean(ensembleMI_control, 2);
sum_wn = nanmean(sum_single_wn, 2);
sum_control = nanmean(sum_single_control,2);
k_wn = sum_wn(end)./n_ensemble(end);
k_control = sum_control(end)./n_ensemble(end);

e_wn = nanstd(ensembleMI_wn,0,2)/sqrt(n_trials);     %Find ste across trials
e_control = nanstd(ensembleMI_control,0,2)/sqrt(n_trials); 
es_wn = nanstd(sum_single_wn,0,2)/sqrt(n_trials);     %Find ste across trials
es_control = nanstd(sum_single_control,0,2)/sqrt(n_trials); 

figure(1);
errorbar(n_ensemble,average_wn,e_wn, 'r');
hold on;
errorbar(n_ensemble, sum_wn, es_wn, 'r--');
plot([0 n_ensemble(end)], [0 sum_wn(end)]);
errorbar(n_ensemble,average_control,e_control,'k');
errorbar(n_ensemble, sum_control, es_control, 'k--');
plot([0 n_ensemble(end)], [0 sum_control(end)], 'k');
legend('WN Ensemble', 'WN Sum', 'Fit', 'Control Ensemble', 'Control Sum', 'Fit');
xlabel('Ensemble size');
ylabel('Information (bits/s)');
hold off;

% Ensemble information per spike
averagespike_wn = nanmean(ensembleMIspike_wn,2);               %Find average information across trials
averagespike_control = nanmean(ensembleMIspike_control, 2);


espike_wn = nanstd(ensembleMIspike_wn,0,2)/sqrt(n_trials);     %Find ste across trials
espike_control = nanstd(ensembleMIspike_control,0,2)/sqrt(n_trials); 


ph = figure(2);
errorbar(n_ensemble,averagespike_control,espike_control,'k','LineWidth', 1);
hold on;
errorbar(n_ensemble,averagespike_wn,espike_wn, 'r', 'LineWidth', 1);
% errorbar(n_ensemble, sum_wn, es_wn, '--');
% plot([0 n_ensemble(end)], [0 sum_wn(end)]);
% errorbar(n_ensemble,average_control,e_control,'r');
% errorbar(n_ensemble, sum_control, es_control, 'r--');
% plot([0 n_ensemble(end)], [0 sum_control(end)], 'r');
% legend('WN Ensemble', 'WN Sum', 'Fit', 'Control Ensemble', 'Control Sum', 'Fit');
xlabel('Ensemble size');
ylabel('Information (bits/spike)');
hold off;

% Redundancy plot
average_wn = nanmean(redundancy_wn,2);
average_control = nanmean(redundancy_control,2);
e_wn = nanstd(redundancy_wn,0,2)/sqrt(n_trials);
e_control = nanstd(redundancy_control,0,2)/sqrt(n_trials);

figure(3);
errorbar(n_ensemble,average_control,e_control,'k','LineWidth',1);
hold on;
errorbar(n_ensemble,average_wn,e_wn, 'r','LineWidth',1);
legend('Control','WN');
xlabel('Ensemble size')
ylabel('Redundancy');
hold off;