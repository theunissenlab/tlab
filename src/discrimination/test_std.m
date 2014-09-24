function [raw_stdSelf raw_stdOther mean_stdSelf mean_stdOther] = test_std(birdname, brainregion, cellname, stimtype)
%[nfiles, stimrate, stdstimrate, backrate, stdbackrate, avgzscore, stdzscore, avgdprime, stddprime, percorrect, avgrank, info, noiseinfo, rateinfo, info1, info2, info2Interval]
% Calculates 3 measures of neural discrimination, the within dprime, the
% percent correct of an ideal oberver based on VR method and the
% information based on Gamma model of spiking neuron
% read all the spike arrival times in a folder

% Read input data from data base
[nfiles spike_times stim_len] = read_all_spikes(birdname, brainregion, cellname, stimtype);

% Window sizes for calculations that depend on window length
winSize = [10 30 50 100];    % Window size for static Info calculation and for ideal observer
ns = length(winSize);

raw_stdSelf = cell(1, ns);
raw_stdOther = cell(1, ns);
mean_stdSelf = cell(1, ns);
mean_stdOther = cell(1, ns);

if nfiles
    
    % Calculate the stds
    for is=1:ns
        [stdSelf stdOther] = calc_std(nfiles, spike_times, stim_len, @VR_distance, winSize(is));

        raw_stdSelf{1, is} = stdSelf;
        raw_stdOther{1, is} = stdOther;
        
        mean_stdSelf{1, is} = mean(stdSelf, 2);
        mean_stdOther{1, is} = mean(stdOther, 2);
        
        fprintf(2, num2str(is));
    end
    
end

return




