function out = faster_surprise_of_one_frequency_band_7D(f_band,stims,stim_names,stim_familiarity,stim_used,w,d,b)
%  Returns the surprise of the stimulus spectrogram band "f_band" of all
%  stimuli in "stims", given the
%  corpus of stims in a cell array with N elements of of B x T spectrograms
%  (B is number of bands, T is time)
%  "stims," their names in cell array "stim_names", their associated
%  familiarity "stim_familiarity" in a 1 x N vector of relative familiarity
%  (use 0 if you don't want a particular stimulus to be included in the
%  corpus), "stim_used" as a row of 0s and 1s corresponding to whether the
%  stimulus will be one presented and a domain described by w (the domain
%  width in ms), d (the minimum gap between stimulus predicted and the domain), and b (the
%  number of stimulus bands to include on either side of the stimulus band
%  predicted).  "out" is a structure of the surprisingly loud and quiet
%  elements associated with the input stims.

%  By default, assume all stims might be used for surprise-STRFing:

if ~exist('stim_used','var')
    stim_used = ones(1,size(stims));
end

%  Default D parameters from Gill et. al. 2008:

if ~exist('w','var')
    w = 3;
end

if ~exist('d','var')
    d = 3;
end

if~exist('b','var')
    b = 4;
end

big_spect = [];
big_familiarity = [];
use_stim = [];
lengths = [];
for istim = 1:length(stims)  %Concatenate the spectrograms
    big_spect = [big_spect stims{istim}];
    big_familiarity = [big_familiarity ones(1,size(stims{istim},2))*stim_familiarity(istim)];
    use_stim = [use_stim ones(1,size(stims{istim},2))*stim_used(istim)];
    lengths = [lengths size(stims{istim},2)];
end
disp('Getting PCA of data...')
tic;
try
    [PCA_Tmat,E,D] = surprise_preprocess(big_spect,f_band,w,d,b); %Transform into the PC space + target of the stimuli
catch
    disp('Problem with spect2PCA_data3:probably out of range.')
    disp(lasterr);
    disp('Exiting rather than erroring...');
    return
end
toc;
PC_Nbins = [15 4 6 6 6 6 15]; % How many bins to use for S and each kept PC.  The first entry is for S, and the next is for the *weakest* kept PC.
n_dims_minus_1 = length(PC_Nbins)-2;

temp_PCA_mat = PCA_Tmat([1 (end-n_dims_minus_1):end],:,1);  % Just the parts of PCA_Tmat that will be used.
clear PCA_Tmat  %Free up some memory
PC_Maxes = max(temp_PCA_mat,[],2);
PC_Mins = min(temp_PCA_mat,[],2);
disp('  ....  finding joint_p P(S,D) table');
tic;
[joint_p,which_cubes,what_position] = surprise_get_joint_p_mat(temp_PCA_mat,PC_Maxes,PC_Mins,PC_Nbins,big_familiarity);
toc;
tic;
disp('  ....  finding cond_p P(S|D) table');
cond_p = surprise_get_7D_cond_p(joint_p);
clear joint_p;
toc;
disp('  ....  finding the probability of each stim')
tic;
out = faster_surprise_cond_p_7D_to_output(cond_p,which_cubes,what_position,use_stim,lengths,stim_names);
toc;
disp('Done this frequency band.');

