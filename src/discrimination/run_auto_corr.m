load 'C:\Documents and Settings\Frederic\My Documents\Matlab\STRF\Analysis\split_dump.mat';


% Auto Correlation functions for field L
stimtype = 'conspecific';
brainregion = 'maybe_L';

n_L = length(L_used_names);

for iL=1:n_L
    [birdname cellname]=strtok(L_used_names{iL},'_');
    cellname = cellname(2:end);
    output_L(iL) = spike_auto_correlation(birdname, brainregion, cellname, stimtype);
end

% Auto Correlation functions for field L
stimtype = 'conspecific';
brainregion = 'maybe_mld';

n_mld = length(mld_used_names);

for imld=1:n_mld
    [birdname cellname]=strtok(mld_used_names{imld},'_');
    cellname = cellname(2:end);
    output_mld(imld) = spike_auto_correlation(birdname, brainregion, cellname, stimtype);
end

clear stimtype birdname brainregion cellname;
save 'C:\Documents and Settings\Frederic\My Documents\Matlab\STRF\Analysis\split_dump_wrate.mat';

% quick examination of isi distribution function
p1_theory_mld = expcdf(1,mean([output_mld.stimrate]));
p100_theory_mld = expcdf(100,mean([output_mld.stimrate]));
p1_100_theory_mld = p1_theory_mld./p100_theory_mld;
fprintf(1,'Actual MLd p1 =%f Poisson = %f\n', mean([output_mld.p1_100ms]), p1_100_theory_mld);

p1_theory_L = expcdf(1,mean([output_L.stimrate]));
p100_theory_L = expcdf(100,mean([output_L.stimrate]));
p1_100_theory_L = p1_theory_L./p100_theory_L;
fprintf(1,'Actual L p1 =%f Poisson = %f\n', mean([output_L.p1_100ms]), p1_100_theory_L);