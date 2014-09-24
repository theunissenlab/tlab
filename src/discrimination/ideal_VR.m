function [percorrect, avgrank] = ideal_VR(nfiles, spike_times, stim_len, winSize)
% Calculates the percent correct detection for one song based on a single trial and the VR metric

nVRruns = 10;       % Number of random choices
debugfigs = 0;      % Set to zero or one for getting debugging figures
template_npts = 1000;           % Use the first 1.0 seconds as a template

ntau = length(winSize);
percorrect = zeros(1, ntau);
avgrank = zeros(1,ntau);

for itau=1:ntau

    tauVR = winSize(itau);         % ms for the exponential smoothing


    %template_npts = ceil(min(stimlen));

    if (template_npts > min(stim_len) )
        fprintf(1, 'Warning in indeal_VR: the template length (%d) is longer than the smallest stimulus (%d ms)\n', ...
            template_npts, min(stim_len));
    end

    expwav = 0:5*round(tauVR);
    expwav = exp(-expwav./tauVR);
    exp_npts = length(expwav);

    % add space for 30 ms offset at begining and end of song
    test_beg = 0;
    test_npts = template_npts+1;
    diff_npts =  test_npts - template_npts;
    correct_ans = 0;
    wrong_ans = 0;
    ntests = 0;
    avgrank(itau) = 0;

    for iruns=1:nVRruns
        % Pick a random spike templates for each stimulus
        template_ind = zeros(1,nfiles);
        template_wav = zeros(nfiles, template_npts+exp_npts-1);
        template_dist = zeros(1, nfiles);
        for nf=1:nfiles
            nt = length(spike_times{nf});
            template_ind(nf) = ceil(nt*rand(1,1));
            time_ind = round(spike_times{nf}{template_ind(nf)});
            nspikes = length(time_ind);
            template_spikes = zeros(1, template_npts);
            for is=1:nspikes
                if (time_ind(is) > 0 && time_ind(is) <= template_npts)
                    template_spikes(1,time_ind(is)) = 1;
                end
            end
            template_wav(nf,:) = conv(template_spikes, expwav);
            template_dist(nf) = sum(template_wav(nf,:).^2);
        end

        if (debugfigs)
            figure(1);
            imagesc(template_wav);
            pause;
        end

        % For all spike trains in nf1 but template find best match by looping through all other.
        for nf1=1:nfiles
            nt = length(spike_times{nf1});
            for it=1:nt
                if it == template_ind(nf1)
                    continue;
                end
                time_ind = round(spike_times{nf1}{it}) + test_beg;
                nspikes = length(time_ind);
                test_spikes = zeros(1, test_npts);
                for is=1:nspikes
                    if (time_ind(is) > 0 && time_ind(is) <= test_npts)
                        test_spikes(1,time_ind(is)) = 1;
                    end
                end
                test_wav = conv(test_spikes, expwav);

                distance_to_template = zeros(1, nfiles);
                offset_to_template = zeros(1, nfiles);

                for nf2=1:nfiles
                    dvals = zeros(1, diff_npts);
                    for dt=1:diff_npts
                        dvals(1,dt) = sum((test_wav(1,dt:template_npts+exp_npts+dt-2)-template_wav(nf2,:)).^2);
                    end
                    distance_to_template(nf2) = min(dvals);
                    offset_to_template(nf2) = find(dvals == distance_to_template(nf2), 1);
                    %distance_to_template(nf2) = (template_dist(nf2)- distance_to_template(nf2))./template_dist(nf2);
                end

                best_distance = min(distance_to_template);
                nf_choice =  find(distance_to_template == best_distance);
                [sorted_distances, sorted_ind] = sort(distance_to_template,'ascend');
                %fprintf(1,'Testing trial %d for stim %d. Best match is stim %d:\n', it, nf1, sorted_ind(1));
                if ( sorted_ind(1) == nf1 )
                    correct_ans = correct_ans+1;
                    rank_correct = 1;
                    % fprintf(1, '\tHit:');
                else
                    wrong_ans = wrong_ans+1;
                    rank_correct = find(sorted_ind == nf1);
                    % fprintf(1, '\tMiss: correct answer has rank %d', rank_correct);
                end
                % fprintf(1, 'Correct %d Wrong %d Percent %f\n', correct_ans, wrong_ans, correct_ans./(correct_ans+wrong_ans));
                avgrank(itau) = avgrank(itau) + rank_correct;
                ntests = ntests + 1;

                % Debugging figure
                if (debugfigs)
                    figure_wav = zeros(nf+1, test_npts+exp_npts-1);
                    figure_wav(1,:) = test_wav;
                    for nf2=1:nfiles
                        figure_wav(1+nf2, offset_to_template(nf2):offset_to_template(nf2)+template_npts+exp_npts-2) = template_wav(nf2,:);
                        figure(nf2+1);
                        plot(figure_wav(1,:));
                        hold on;
                        plot(figure_wav(nf2+1,:),'r');
                        title(sprintf('%f', distance_to_template(nf2)));
                        hold off;
                    end
                    pause;
                end
            end
        end
    end

    percorrect(itau) = correct_ans./(correct_ans+wrong_ans);
    avgrank(itau) = avgrank(itau)./ntests;

end

return