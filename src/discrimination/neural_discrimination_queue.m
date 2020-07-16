function neural_discrimination_queue(birdname, brainregion, cellname, stimtype)
% Runs neural discrimination and saves the results in a mat file

outputs = neural_discrimination(birdname, brainregion, cellname, stimtype);
% save to a file
done_name = sprintf('%s_%s_%s_%s.mat' , birdname, brainregion, cellname, stimtype);
save(done_name);

return




