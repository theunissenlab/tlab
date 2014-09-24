% birdname = 'gg0304';   % This is a very good neuron for song discrimination
% cellname = '6_B'; 
birdname = 'blahpi0303';  % This is an average neuron
cellname = '10_B';
stimtype = 'conspecific';
brainregion = 'maybe_L';

%output = spike_auto_correlation(birdname, brainregion, cellname, stimtype);

outputs4 = neural_discrimination_limited(birdname, brainregion, cellname, stimtype, 4, 10);
outputs8 = neural_discrimination_limited(birdname, brainregion, cellname, stimtype, 8, 10);
outputs12 = neural_discrimination_limited(birdname, brainregion, cellname, stimtype, 12, 10);
outputs16 = neural_discrimination_limited(birdname, brainregion, cellname, stimtype, 16, 10);
outputs20 = neural_discrimination_limited(birdname, brainregion, cellname, stimtype, 20, 1);

save zfnormal_fieldL_blahpi0303_10_B.mat outputs4 outputs8 outputs12 outputs16;