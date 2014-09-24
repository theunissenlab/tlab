function new_song = songfilt_call(song, fs, f_low, f_high, db_att, call_type)
% Filters and normalizes song between -1 and 1 based on mean power in non silence regions.
% Use db_att to get smaller values of power


% Deal with arguments and parameters
if isempty(song)
    error('Missing or empty song array'); 
end
if isempty(fs)
    error('Missing or empty song array'); 
end

if isempty(f_low)
    f_low = 250.0;
end

if isempty(f_high)
    f_high = 8000.0;
end

NUMBER_STD = 3;
 

if isempty(db_att)
    db_att = 0;
end

if  db_att == 0.0
    extmult = 1.0;
else 
    extmult = 10.0.^(db_att/-20.0);
end

% set up Lindex the index that make power attenuation to respect
% differences of loudness between call categories. Values calculated from
% the mean values in the data bank of first recording group. calculated
% using cal_Mean_power on each series of calls (1 series of calls per
% category)
Lindex = 1;
if strcmp(call_type,'Ag')
    Lindex = 0.0297/0.186;
end

if strcmp(call_type,'Be')
    Lindex = 0.0228/0.186;
end

if strcmp(call_type,'DC')
    Lindex = 0.0854/0.186;
end

if strcmp(call_type,'DCST')%power of non-stressed and stressed DC from saint-etienne
    Lindex = 0.0714/0.186;
end

if strcmp(call_type,'Di')
    Lindex = 0.0495/0.186;
end

if strcmp(call_type,'LT')
    Lindex = 0.0172/0.186;
end

if strcmp(call_type,'Ne')
    Lindex = 0.0160/0.186;
end

if strcmp(call_type,'So')
    Lindex = 1;
end

if strcmp(call_type,'Te')
    Lindex = 0.0100/0.1806; 
end

if strcmp(call_type,'Th')
    Lindex = 0.0333/0.1806;
end

if strcmp(call_type,'Wh')
    Lindex = 0.0069/0.1806;
end


% Figure out length of filter
nframes = length(song);
if ( nframes > 3*512 ) 
    nfilt = 512;
elseif ( nframes > 3*64 ) 
    nfilt = 64;
    
elseif ( nframes > 3*16 ) 
    nfilt = 16;
else
    error('Song data is too short for filtering');
end

% Generate filter and filter song
song_filter = fir1(nfilt,[f_low*2.0/fs, f_high*2.0/fs]);
new_song = filtfilt(song_filter, 1, song);

% Rescale file to get maximum dynamic range 

max_song = max(new_song);
min_song = min(new_song);
range = max(-min_song,max_song)/10.0;

pow =0.0;
nclip = 0;
for i = 2:nframes-1 
    val = new_song(i);
    if ( val > range & val > new_song(i-1) & val > new_song(i+1) )
        nclip=nclip+1;
        pow = pow + val*val;
    end
    if ( val < -range & val < new_song(i-1) & val < new_song(i+1) )
        nclip=nclip+1;
        pow = pow+val*val;
    end
end
pow = sqrt(pow/nclip);
range = pow*NUMBER_STD;

new_song = new_song .* (extmult*Lindex);
nclip = 0;
for i =1:nframes 
    
    if ( new_song(i) > range ) 
        new_song(i) = range;
        nclip=nclip+1;
        
    elseif ( new_song(i) < -range )
        new_song(i) = -range;
        nclip=nclip+1;
    end
end
new_song = new_song./range;

if ( nclip ~= 0 )
    warning_str = sprintf('Warning: %d points were clipped', nclip);
    warning(warning_str);
end

% End of function songfilt_call.m
