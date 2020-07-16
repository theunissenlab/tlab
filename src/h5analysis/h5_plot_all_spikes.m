function h5_plot_all_spikes(h5Path, stimType, pThreshold, spikeSnipFlg, spikeSnipFile)
% Plots spectrogram and spike rasters for all stimuli of one type recorded
% at one site.
% stimType is the class of stimulus played - usually a protocol name such
% as Mask1, Call1, etc.
% pThreshold is a floating point value to plot only responses with a p
% value below a cuttof.  To plot all set pThreshold to 1.0.
% spikeSnipFlg is using to also plot the snippets for all spikes.

unit = read_unit_h5file(h5Path, 'r');

% If plottin spike shapes get pointer to snippets
if exist('spikeSnipFlg')
    % Open Snipet file
    C=textscan(h5Path, '%s', 'delimiter','/');
    blockname = C{1}{end};    % we might be able to use this at the end.
    
    fid = fopen(spikeSnipFile, 'r');
    if fid == -1
        fprintf(1, 'Error opening waveform snips files %s\n', spikeSnipFile);
        return;
    end
else
    fid = -1;
end

nclasses = length(unit.classes);

classId = 0;
for ic=1:nclasses
    if strcmp(stimType, unit.classes{ic})
        classId = ic;
        break;
    end
end
if classId == 0
    fprintf(1, 'Warning: could not find stimType %s in h5 file %s\n', stimType, h5Path);
    fprintf(1, '\tAvailable options are:\n');
    for ic=1:nclasses
        fprintf(1,'\t\t%s\n', unit.classes{ic});
    end
    return;
end

responses = unit.class_responses.(stimType);

% This is the number of sound files played
nfiles = length(responses);


ngood = 0;
for nfi=1:nfiles
    if responses{nfi}.pvalue < pThreshold
       h5_plot_one_spike(responses{nfi}, fid);
       ngood = ngood + 1;
       pause();
    end
end
fprintf(1,'Found %d stim out of %d with significant responses (p<%.3f)\n', ngood, nfiles, pThreshold);


    



return




