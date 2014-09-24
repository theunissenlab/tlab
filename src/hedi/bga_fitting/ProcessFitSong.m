function fsong=ProcessFitSong(r_song)
% process whole fitted song after being fitted
% 

% first compute the whole waves 

fsong=SPLtoWave(r_song);

% smooth the song 
%
%fsong=SmoothSong(fsong);

% compute waves for smoothed 

%fsong=SmoothSPLtoWave(fsong);

% compute plot and display and whole length motograms

fsong=SPLtoWholeMotograms(fsong);

