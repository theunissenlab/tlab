function check_signif_spikes(h5Path, stimType, pThreshold)
% Gives out the number of significant stims over the total number of stims
% for the whole h5 file, depending on the chosen pThreshold.
% stimType is the class of stimulus played - usually a protocol name such as Mask1, Call1, etc.
% pThreshold is a floating point value to plot only responses with a p value below a cuttof.

unit = read_unit_h5file(h5Path);
nclasses = length(unit.classes); % number of protocols (ex. Mask1)
responses = unit.class_responses.(stimType);
nresp = size(responses,2);

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

% This is the number of sound files played
nfiles = length(responses);

ngood = 0;
for nfi=1:nfiles
    if responses{nfi}.pvalue < pThreshold
       ngood = ngood + 1;
    end
end
fprintf(1,'Found %d stim out of %d with significant responses (p<%.3f)\n', ngood, nfiles, pThreshold);

return
