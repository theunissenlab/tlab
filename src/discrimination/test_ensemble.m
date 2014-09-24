ncells = 1;

% cellid(1).birdname = 'pipu1018';
% cellid(1).cellname = '4_B'; 
% cellid(1).brainregion = 'maybe_L';
cellid(1).birdname = 'pipu1018';
cellid(1).cellname = '6_B'; 
cellid(1).brainregion = 'maybe_L';

stimtype = 'conspecific';
output_ensemble = neural_discrimination_ensemble(cellid, stimtype, 1000, []);
