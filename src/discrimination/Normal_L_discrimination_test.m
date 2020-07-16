% Runs the neural discrimination for normal L
brainregion = 'maybe_L';
fname = 'zfnormal_fieldL.txt';

% Load all cells that share 20 songs in field L
load birds.mat

% Allocate space for output values
number_cells = length(birds);

% Check to see if all files exists...
% First conspecific
stimtype = 'conspecific';

% for testing first
number_cells = 1;

% for nc=1:number_cells
nc = 10;
    [birdname remain]= strtok(birds{nc}, '_');
    [cellname remain]= strtok(remain, '_');
    cellname = strcat(cellname, remain);
    outputs = neural_discrimination(birdname, brainregion, cellname, stimtype);
% end

