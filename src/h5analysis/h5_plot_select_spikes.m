function h5_plot_select_spikes(h5Path, stimType, pThreshold, inputarg)
% Plots spectrogram and spike rasters for selected stimuli of one type recorded at one site.
% inputarg is a cell array of cells with the selected parameters, eg. {{'birdid', 'f6'}, {etc}}.
% stimType is the class of stimulus played - usually a protocol name such as Mask1, Call1, etc.
% pThreshold is a floating point value to plot only responses with a p value below a cuttof.  
% To plot all set pThreshold to 1.0.
% For now, only one value can be entered for each argument (eg. 'f4').

% ex of h5Path: h5Path = '/auto/fdata/solveig/ELECTROPHY/h5/BlaLb1986M/Site1_R2408_7.h5';
% ex of inputarg : inputarg = {{'birdid', 'f15'}, {'callnum', '12'}};

% check presence of arguments
nargs = size(inputarg,2);

if nargs == 0
    fprintf('You did not select any selective argument for this function. Please use the script h5_plot_all_spikes.m instead.\n');
return;
end

% retrieve data
unit = read_unit_h5file(h5Path, 'r');
nclasses = length(unit.classes); % number of protocols (ex. Mask1)
fid = -1; % for now this script doesn't plot snipets (see line 82)    !! A VERIFIER !!!

classId = 0;
for ic=1:nclasses
    if strcmp(stimType, unit.classes{ic})
        classId = ic;
        break;
    end
end

if classId == 0  % Check for the available parameters
    fprintf(1, 'Warning: could not find stimType %s in h5 file %s\n', stimType, h5Path);
    fprintf(1, '\tAvailable options are:\n');
    for ic=1:nclasses
        fprintf(1,'\t\t%s\n', unit.classes{ic});
    end
    return;
end

responses = unit.class_responses.(stimType);
nresp = size(responses,2);
fields = fieldnames(responses{1}); % available parameters

% further check the arguments
if ~iscell(inputarg{1}) %only on eargument!
    nargs=1;
    argName={inputarg{1}};
    argVal={inputarg{2}};
else
    argName = {}; argVal = {};
    for i = 1:nargs
        if size(inputarg{i}, 2) ~= 2
            fprintf(2,'Invalid inputarg #%d. Enter a vector of cells, each cell containing\n', i);
            fprintf(2, 'the argument name (eg. ''birdid'') and the argument value (eg. ''f4'').\n');
            return;
        elseif isempty(find(strcmp(fields, inputarg{i}(1)), 1))
            fprintf(2, 'Invalid inputarg name #%d.\n', i);
            fprintf(1, 'Available options are:\n');
            disp(fields)
            return;
        else argName = [argName inputarg{i}(1)];
            argVal = [argVal inputarg{i}(2)];
        end
    end
end

% loop through the arguments and select the corresponding stimuli
stimSel = []; % this vector will contain the numbers of the selected stimuli on the TDT
stimNum = [];
for j = 1 : nresp
    fieldVal = {};
    for k = 1 : nargs
        Val = getfield(responses{j}, argName{k});
        fieldVal = [fieldVal Val];
    end
    Cmp = strcmp(fieldVal, argVal);
    check = find(Cmp == 0);
    if isempty(check)
        ind = getfield(responses{j}, 'number');
        stimNum = [stimNum str2num(ind)]; % This contains the actual stim numbers
        stimSel = [stimSel j];
    end
end

nfiles = length(stimSel); % this is the number of sound files selected
if nfiles == 0
    fprintf('No stimulus corresponds to the parameters you selected.\n');
    return;
end

% plot selected files
ngood = 0;
for nfi = stimSel
    if responses{nfi}.pvalue < pThreshold
       figure()
       h5_plot_oneSel_spike(responses{nfi}, fid);
       ngood = ngood + 1;
       %pause();
    end
end
fprintf(1,'Found %d stim with significant responses out of the %d selected (p<%.3f)\n', ngood, nfiles, pThreshold);

return


