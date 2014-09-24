[song, samprate]= wavread('/Users/frederictheunissen/Documents/Data/birdsongs/ucsf/zfa_25/zfa_25.6.filt.wav'); 

% These are the parameters we used to get our STRFs
nstd = 6;
fband = 125;
twindow = nstd/(fband*2.0*pi);   % Window length
winLength = fix(twindow*samprate);  % Window length in number of points
winLength = fix(winLength/2)*2;


flow = 500;
fhigh = 8000;
dBScale = 50;

figure(1);
spec(song, fband, samprate, dBScale);

[s, fc] = GaussianFilterBank(song, winLength, samprate, flow, fhigh);

for i=1:length(fc)
    if (i < 5)
        figure(i+1);
        spec(s(i,:), fband, samprate, dBScale);
        title(sprintf('Narrow Band signal centered at %f', fc(i)));
    end
end

% Reconstruct sound
songRecon = sum(s);

soundsc(song, samprate);
pause();

soundsc(songRecon, samprate);
figure(6);
spec(songRecon, fband, samprate, dBScale);




