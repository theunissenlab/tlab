function [freq intensity rateDiff dur] = rateDiff_noopur_spikes(group, birdname, cellname)
% Plots spectrogram and spike rasters for all stimuli of one type recorded
% at one site.

nw = 512;   % window for psd and minimum duration stimulus
DELAY = 0.015; % Allow for a 15 ms delay
ANALTIME = 0.050;  % Look at 50 ms after beginning only

[nfiles, spike_times, stimlen, stimname] = read_noopur_spikes(group, birdname, cellname);

freq = [];
intensity = [];
rateDiff = [];
dur = [];

if nfiles
    
    % Get the average background rate
    [stimrate stdstimrate backrate stdbackrate avgzscore stdzscore avgdprime stddprime] = dprime_within(nfiles, spike_times, stimlen);
    
    if (avgzscore - 2*stdzscore/sqrt(nfiles) > 0)
        for nfi=1:nfiles
            spike_time = spike_times{nfi};
            stim_name = stimname{nfi};
            
            % Get sound, number of trials, etc
            ntrials = length(spike_time);
            
            [sound_in, samprate , nbits] = wavread(stim_name);
            
            
            % Divide the sound into areas of non-zero time
            ind = find(sound_in);
            on = zeros(size(sound_in));
            sound_p1 = zeros(size(sound_in));
            sound_p1(1:end-1) = sound_in(2:end);
            sound_m1 = zeros(size(sound_in));
            sound_m1(2:end) = sound_in(1:end-1);
            sound_diff = sound_p1-sound_m1;
            ind2 = find(sound_diff);
            
            % Sounds is "on" when value is not zero or derivative is not zero
            on(ind) = 0.8;    % It is set at 0.8 to make a nice graphic.
            on(ind2) = 0.8;
            
%             plot(sound_in);
%             hold on;
%             plot(on, 'r');
%             hold off;
%             pause();
            
            % Find beggining and end
            laston = 0;
            for i=1:length(on)
                if (laston == 0) && (on(i)~=0)
                    begInd = i;
                elseif (laston~=0) && (on(i) == 0)
                    endInd = i;
                    durInd = endInd - begInd;
                    
                    if (durInd > nw )    % Only examine stimulus intervals of certain length
                        % calculate duration, peak frequency and intensity
                        durVal = (endInd-begInd)/samprate;
                        intVal = std(sound_in(begInd:endInd));
                        [Pxx,f] = pwelch(sound_in(begInd:endInd), nw, fix(nw/2) , nw, samprate);
                        Pmax = max(Pxx);
                        fVal = f(Pxx==Pmax);
                        
                        % calculate number of spikes during interval
                        nspikes = 0;
                        begTime = begInd/samprate;
                        endTime = endInd/samprate;
                        
                        for it=1:ntrials
                            spike_time_trials = spike_time{it}./1000;    % Spike time are in ms and code is in s
                            ns = length(spike_time_trials);                            
                            for is=1:ns
                                if (spike_time_trials(is) >= begTime)
                                    %if (spike_time_trials(is) <= (endTime + DELAY) )
                                    if (spike_time_trials(is) <= (begTime + ANALTIME) )
                                        nspikes = nspikes + 1;
                                    else
                                        break;
                                    end
                                end
                            end
                        end
                        % stimSectionRate = nspikes/(durVal*ntrials);
                        stimSectionRate = nspikes/(ANALTIME*ntrials);
                        rateDiffVal = stimSectionRate - backrate;
                        

                        % Stuff values
                        intensity = [intensity intVal];
                        freq = [freq, fVal];
                        rateDiff = [rateDiff rateDiffVal];
                        dur = [dur durVal];
                        % fprintf(1, 'Beg= %.2f End= %.2f Int= %.2f Dur= %.2f Freq= %.0f RD= %.2f\n', begTime, endTime, 20*log10(intVal), durVal, fVal, rateDiffVal);
                    end
                end
                laston = on(i);
            end
            
            
        end
    end
end

return




