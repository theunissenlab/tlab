function [rasters, FrameRate_Hz] = tmpraster(PF)
%function [rasters, FrameRate_Hz] = plxraster(fname, stimsettings)
%
% Make binned rasters from a pype-plexon combined data
% The pype data should be compatible with sglshow task.
% (e.g., fsglshow)
%
% See also: addplxspikes2pype.m
%

%disp(['File: ' fname]);

SettleDuration_ms = 200;

%load(fname, 'PF');
FrameRate_Hz = PF.rec(1).params.mon_fps;
frame_ms = 1/FrameRate_Hz*1000;
numtrials = size(PF.rec, 2);

SettleFrames = ceil(SettleDuration_ms/frame_ms);

MaxFrame = NaN;
%% get movie information
for ii=1:size(PF.extradata,2);
    a = cell2mat(PF.extradata(ii));
    if strmatch(a.id, 'header movie')
        b = strvcat(a.data{1});
        frame_i = strmatch('# framecount', b);
        [dummy  dummy MaxFrame] = strread(b(frame_i,:), '%s%s%d', 'delimiter', ' ');
    end
end

if isnan(MaxFrame)
    error('Can not find frame number information in the p2m file!');
end

% fprintf('Number of frames: %d\n', MaxFrame);
% fprintf('Frame rate : %d Hz\n', FrameRate_Hz);
% fprintf('Number of trials : %d Hz\n', numtrials);


nu = []; nc = [];
for ii = 1:numtrials
    nu= [nu size(PF.rec(ii).plexon_spikes_final,1)];
    nc= [nc size(PF.rec(ii).plexon_spikes_final,2)];
end

nu = max(nu); nc = max(nc);

rasters=cell(nu,nc);

rasters(:) = {nan([MaxFrame numtrials], 'single')};

for tt = 1:numtrials
    if isempty(PF.rec(tt).plexon_spikes_final)
        PF.rec(tt).plexon_spikes_final = cell(nu,nc);
    end
    % for ui=1:nu
    %     for ci=1:nc
    %         if ~isempty(PF.rec(tt).plexon_spikes) || ~isempty(PF.rec(tt).plexon_spikes{ui,ci})
    %             % rasters{ui,ci} = single(repmat(NaN, [MaxFrame numtrials]));
    %             rasters{ui,ci} = nan([MaxFrame numtrials], 'single');
    %         end
    %     end
    % end
end

% figure(1); clf;
for trialnum=1:numtrials
    displayed = 0;

    thisrec = PF.rec(trialnum);

    if isempty(thisrec.plexon_spikes_final)
        continue;
    end
    % obtain flip information
    flip_i = strmatch('frame_flipped', thisrec.ev_e);
    if length(flip_i)<10, continue, end;

    flip_t = thisrec.ev_t(flip_i);

    for fi=1:length(flip_i)
        [dummy, fnum, pol] = strread(cell2mat(thisrec.ev_e(flip_i(fi))), '%s%d%s', 'delimiter',':');
        flip_num(fi) = fnum;
        if strmatch(pol, 'True')
            flip_polarity(fi) = 1;
        elseif strmatch(pol, 'False')
            flip_polarity(fi) = 0;
        else
            flip_polarity(fi) = NaN;
        end
    end

    % obtain blit information
    blit_i = strmatch('image blit', thisrec.ev_e);
    blit_t = thisrec.ev_t(blit_i);

    for bi=1:length(blit_i)
        [dummy, fnum] = strread(cell2mat(thisrec.ev_e(blit_i(bi))), '%s%d', 'delimiter',':');
        blit_num(bi) = fnum;
    end

    % obtain task timings
    evs = thisrec.ev_e;
    trial_done_i = [0 strmatch('fix_done', evs) strmatch('bar_up', evs) strmatch('fix_lost', evs) strmatch('eye_stop', evs)];
    trial_done_i = trial_done_i(2:end);
    possibledonetime = [thisrec.ev_t([trial_done_i]) flip_t(end)];
    trial_done_time = min(possibledonetime);

    % obtain photo timings
    st = thisrec.ev_t(flip_i(1));
    photo_time = thisrec.photo_t;

    % obtain frame drops
    photo_bin = 0:frame_ms:photo_time(end)-photo_time(1)+frame_ms;
    fp=hist(photo_time-photo_time(1), photo_bin);

    fstep = ones(size(fp));
    for fpi=2:length(fp)
        if fp(fpi) == fp(fpi-1)
            fstep(fpi)=0;
        end
    end
    fdrops = (fstep==0);
    fcum = cumsum(fstep);
    first_frame = flip_num(2);
    fframe = fcum+first_frame-1;
    fframe = mod(fframe-1,MaxFrame)+1;

    % mask initials and frame drops
    masksize = round(SettleDuration_ms/frame_ms);
    fframe(1:masksize) = NaN;
    fmax = length(fcum);
    for ii=find(fdrops)
        masks = ii:masksize;
        masks = find(masks<fmax);
        fframe(masks) = NaN;
    end
    if sum(fdrops)
        fprintf('d%d',sum(fdrops));
    else
        fprintf('o');
    end
    validbins = ~isnan(fframe);

    % plot some info
%     subplot(6,1,mod(trialnum-1,6)+1);
%     cla;
%     if isfield(thisrec, 'plexon_photos')
%         plot(diff(thisrec.plexon_photos),'rx-');
%     end
%     hold on;
%     plot(diff(photo_time),'bo-');
%     ylim([10 50]);
%     title(sprintf('Trial %d, dropped %d/%d frames', trialnum, sum(fdrops),length(fframe)));
    
    rtoffset = thisrec.realt(1);
    
    nc1 = size(thisrec.plexon_spikes_final,2);
    
    % bin spikes

        for ci=1:nc1
            if ~isempty(thisrec.plexon_spikes_final) || ~isempty(thisrec.plexon_spikes_final{ci})
                thisspikes = thisrec.plexon_spikes_final{ci};
                thisspikes = thisspikes-rtoffset-photo_time(1);
                thisspikes = thisspikes(thisspikes>=0 & thisspikes<=photo_bin(end));
                sp=hist(thisspikes, photo_bin);
                rasters{ci}(fframe(validbins),trialnum) = sp(validbins);
            end
        end


end

fprintf('\n');

% figure(1);clf;
% for ii=1:32
%     subplot(16,2,ii);
%     imagesc(compress_raster(rasters{1,ii})');
% end
