% Code than randomizes the order of syllables in a song
% Written May 2004 by Frederic Theunissen for Tiny Boumans

% User parameters
songname = 'chain6.wav';
newsongname = 'chain6_rand.wav';
lowf_sound = 1000;              % High pass filtering cutoff for sound
lowf_amp = 100;                % Low pass filtering cuttoff for amplitude enveloppe
threshold = 0.01;

% Read song
[song fs nbits] = wavread(songname);
fnyquist = fs/2;

% High pass filter the song
fnorm = lowf_sound/fnyquist;
b_song = fir1(128, fnorm, 'high');
songfilt = filtfilt(b_song, 1, song);

% Rectify filtered song
songabs = abs(songfilt);


% Low-pass filter to get amplitude enveloppe
fnorm = lowf_amp/fnyquist;
b_amp = fir1(128, fnorm);
songenv = filtfilt(b_amp, 1, songabs);
maxsong = max(songenv);
thval = maxsong*threshold;

plot(song);
hold on;
plot(songenv,'k');
hold on;
plot([1 length(song)], [thval thval], 'k--');

% Find begining and end of syllables
j_beg = 0;
j_end = 0;
for i=1:length(songenv)-1
    if (songenv(i) < thval & songenv(i+1) >= thval )
        j_beg = j_beg + 1;
        beg_ind(j_beg) = i;
        
    elseif (songenv(i) >= thval & songenv(i+1) < thval )
        j_end = j_end + 1;
        end_ind(j_end) = i+1;
    end
end

% Some error checking to make sure that the signal starts with a begining.
if ( beg_ind(1) > end_ind(1) ) 
    beg_ind_temp = beg_ind(2:end);
    clear beg_ind;
    beg_ind = beg_ind_temp;
    clear beg_ind_temp;
end


% Make syllabel boundaries by taking the middle of the silence section.
for j=1:j_beg
    if j == 1
        syll_bound(1,j) = 1;
    else
        syll_bound(1,j) = end_ind(j-1)+ floor((beg_ind(j)-end_ind(j-1))/2);
    end
    if j == j_beg
        syll_bound(2,j) = length(songenv);
    else
        syll_bound(2,j) = end_ind(j) + floor((beg_ind(j+1) - end_ind(j))/2)-1;
    end
    plot([syll_bound(1,j) syll_bound(1,j)], [-maxsong maxsong], 'r');
end

% Reshuffle j_beg index and make new song with syllables in random order
rand_j = randperm(j_beg);
newsong = [];
for j=1:j_beg
    newsong = [newsong; song(syll_bound(1,rand_j(j)):syll_bound(2,rand_j(j)))];
end

% plot new song on new graph
figure;
plot(newsong);

% Write out song 
wavwrite(newsong, fs, newsongname);

