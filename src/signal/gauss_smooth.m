function smoothed = gauss_smooth(input,npoints,nstd)

if nargin < 3
    nstd = 6;
end

win = gausswin(npoints,sqrt(nstd));
win = win./sum(win);

smoothed = conv(input,win,'same');