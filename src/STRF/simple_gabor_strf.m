function strf = simple_gabor_strf(numChannels, delays, x0, y0, params)

    if nargin < 5       
        params = zeros(1, 5);
        params(1) = 0;  %theta
        params(2) = 1;  %temporal bandwidth
        params(3) = .1; %aspect ratio
        params(4) = .1; %frequency bandwidth
        params(5) = 0;  %phase
    end

    numDelays = length(delays);
    
    [xvals, yvals] = meshgrid(0:(numDelays-1), 0:(numChannels-1));    
    strf1 = gabor2d(xvals, yvals, x0, y0, params);
    strf2 = gabor2d(xvals, yvals, x0+3, y0, params);
    strf = strf1 - 1e-1*strf2;
    