% set up path
addpath(genpath('/auto/k2/share/matlab/fet/Neural Discrimination'));
addpath(genpath('/auto/k2/share/matlab/fet/P Data Base Tools'));
addpath(genpath('/auto/k2/share/matlab/fet/Sound Tools'));

% Runs the neural discrimination for normal L
brainregion = 'wnL';
stimtype = 'conspecific';

% These are good
  

% {'pipu1018', '11_B' }, This one is excluded because it is bad

% All the data files pasted from the Excel file
wn_allunits = { ...
    {'wnmale1', '5' }, ... %wnmale1 is 90 days old
    {'wnmale1', '6' }, ... %wnmale1 is 90 days old
    {'wnmale4', '9' }, ... 
    {'wnmale6', '8' }, ... 
    {'wnmale6', '9' }, ... 
    {'wnmale6', '10' }, ... 
    {'wnmale6', '11' }, ... 
    {'wnmale7', '11' }, ... 
    {'wnmale7', '15' }, ... 
    {'wnmale8', '1' }, ... 
    {'wnmale8', '2' }, ... 
    {'wnmale8', '3' }, ... 
    {'wnfemale1', '6' }, ... 
    {'wnfemale1', '7' }, ... 
    {'wnfemale1', '8' }, ... 
    {'wnfemale1', '9' }, ... 
    {'wnfemale4', '12' }, ... 
    {'wnfemale4', '13' }, ... 
    {'wnfemale5', '5' }, ... 
    {'wnfemale5', '7' }, ... 
    {'wnfemale5', '8' }, ... 
    {'wnfemale5', '9' }, ... 
    {'wnfemale5', '12' }, ... 
    {'wnfemale6', '5' }, ... 
    {'wnfemale6', '6' }, ... 
    {'wnfemale7', '13' }, ... 
    {'wnfemale8', '6' }};

number_cells = length(wn_allunits);


for i=1:number_cells
    birdname = wn_allunits{i}{1};
    cellname = wn_allunits{i}{2};
    plot_all_spikes(birdname, brainregion, cellname, stimtype);
    pause();
    close all;
end


