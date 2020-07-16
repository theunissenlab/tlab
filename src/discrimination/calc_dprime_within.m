
birdname = 'gg0304';
cellname = '4_A';
stimtype = 'zfsong';

rootcd = 'C:\Documents and Settings\Frederic\My Documents\Data\STRFs\mld-l-strf\';
savecd = pwd;
datacd = strcat(rootcd, birdname, '\', cellname, '\', stimtype);
cd(datacd);

% Read the list of data files in the folder and the stim length
[datname stimlen]=textread('stim_init.rec','%s %f');
nfiles = length(datname);

% Read all the spike arrival times for each stim 
spike_times = cell(1, nfiles);

for nf=1:nfiles
    % Read the spike arrival data
    remain = datname{nf};
    while true
        [token, remain] = strtok(remain, '/');
        if isempty(remain), break; end
    end
    
    fname = sprintf('%s.txt', token);

    fid = fopen(fname);
    nt = 0;
    while true
        tline = fgetl(fid);
        if ~ischar(tline), break; end
        nt = nt + 1;
    end

    spike_times{nf} = cell(1,nt);
    frewind(fid);
    
    for it=1:nt
        tline = fgetl(fid);
        spike_times{nf}(it) = textscan(tline, '%f');
    end
    
    fclose(fid);

end

% Calculate the average and std spike rate during the stimulus period only.
spike_mean_rate = zeros(1,nf);
spike_std_rate = zeros(1,nf);
for nf=1:nfiles
    nt = length(spike_times{nf});
    spike_count = zeros(1,nt);
    for it=1:nt
        spike_count(it) = length(find(spike_times{nf}{it} > 0.0 & spike_times{nf}{it} < stimlen(nf)));
    end
    spike_mean_rate(nf) = 1000.0*mean(spike_count)./stimlen(nf);
    spike_std_rate(nf) = 1000.0*std(spike_count)./stimlen(nf);
end

% Calculate all the within-dprimes
icomp = 1;
for nf1=1:nfiles-1
    for nf2=nf1+1:nfiles
        dprime(icomp) = 2*abs(spike_mean_rate(nf1)-spike_mean_rate(nf2))./sqrt(spike_std_rate(nf1)^2+spike_std_rate(nf2)^2);
        icomp = icomp +1;
    end
end

% Print out the average d prime
fprintf(1, '%s %s %s %f %f\n', birdname, cellname, stimtype, mean(dprime), std(dprime));

% Calculate VR metric
nVRruns = 10;
tauVR = 30;         % ms for the exponential smoothing
%template_npts = ceil(min(stimlen));
template_npts = 1500;           % Use the first 1.5 seconds.
expwav = 0:5*round(tauVR);
expwav = exp(-expwav./tauVR);
exp_npts = length(expwav);

% add space for 2 second offset and end of song
test_beg = 2000;
test_npts = 6000;

diff_npts =  test_npts - template_npts;

for iruns=1:nVRruns
    % Pick templates
    template_ind = zeros(1,nf);
    template_wav = zeros(nf, template_npts+exp_npts-1);
    template_dist = zeros(1, nf);
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
    figure(1);
    imagesc(template_wav);
    pause();
    
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
                offset_to_template(nf2) = min (find(dvals == distance_to_template(nf2) ));
                distance_to_template(nf2) = (template_dist(nf2)- distance_to_template(nf2))./template_dist(nf2);
            end
            % Debugging figure
            
            figure_wav = zeros(nf+1, test_npts+exp_npts-1);
            figure_wav(1,:) = test_wav;
            for nf2=1:nfiles
                figure_wav(1+nf2, offset_to_template(nf2):offset_to_template(nf2)+template_npts+exp_npts-2) = template_wav(nf2,:);
                figure(nf2+1);
                plot(figure_wav(1,:));
                hold on;
                plot(figure_wav(nf2+1,:),'r');
                title(sprintf('%f', distance_to_template(nf2)));
            end
            pause();
        end
    end
end

cd(savecd);
