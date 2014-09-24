function [stimrate, stdstimrate, backrate, stdbackrate, avgzscore, stdzscore, avgdprime, stddprime, dprime_mat] = dprime_within(nfiles, spike_times, stimlen)
% Calculates the within dprime values by comparing normalized responses for
% all pair-wise combinations.

pre_time = 2000;  % Background time before spike hard coded at 2 s but this should be checked.
% Calculate the average and std spike rate during the stimulus period only
% and the z score obtained by subtracting bacground rate
stim_mean_rate = zeros(1,nfiles);
back_mean_rate = zeros(1,nfiles);
spike_std_rate = zeros(1,nfiles);
spike_zscore = zeros(1,nfiles);
for nf=1:nfiles
    nt = length(spike_times{nf});
    spike_count = zeros(1,nt);
    spike_count_bg = zeros(1,nt);
    for it=1:nt
        spike_count(it) = length(find(spike_times{nf}{it} > 0.0 & spike_times{nf}{it} < stimlen(nf)));
        spike_count_bg(it) = length(find(spike_times{nf}{it} <= 0.0 ));
    end
    stim_mean_rate(nf) = 1000.0*mean(spike_count)./stimlen(nf);
    back_mean_rate(nf) = 1000.0*mean(spike_count_bg)./pre_time;
    spike_std_rate(nf) = 1000.0*std(spike_count)./stimlen(nf);
    drate = spike_count./stimlen(nf) - spike_count_bg./pre_time;
    mdrate = mean(drate);
    if (mdrate == 0 )
        spike_zscore(nf) = 0;
    else
        spike_zscore(nf) = mdrate./std(drate);
    end
end
stimrate = mean(stim_mean_rate);
stdstimrate = std(stim_mean_rate);
backrate = mean(back_mean_rate);
stdbackrate = std(back_mean_rate);
avgzscore = mean(spike_zscore);
stdzscore = std(spike_zscore);

% Calculate all the within-dprimes
icomp = 1;
dprime_mat = zeros(nfiles,nfiles);

for nf1=1:nfiles-1
    for nf2=nf1+1:nfiles
        if ( spike_std_rate(nf1) == 0 && spike_std_rate(nf2) == 0 )
            dprime(icomp) = 0;
        else
            dprime(icomp) = 2*abs(stim_mean_rate(nf1)-stim_mean_rate(nf2))./sqrt(spike_std_rate(nf1)^2+spike_std_rate(nf2)^2);
        end
        dprime_mat(nf1, nf2) = dprime(icomp);
        dprime_mat(nf2, nf1) = dprime(icomp); % attention la diagonale aura des zéros (changer en NA pour calculer la moyenne ?)
        icomp = icomp +1;
    end
end

avgdprime = mean(dprime);
stddprime = std(dprime);


return