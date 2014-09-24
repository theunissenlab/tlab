function val = gabor2d(x, y, x0, y0, params)

    theta = params(1);
    tempBw = params(2);
    aspectRatio = params(3);
    freqBw = params(4);
    phase = params(5);
    
    xp = (x-x0)*cos(theta) + (y-y0)*sin(theta);
    yp = -(x-x0)*sin(theta) + (y-y0)*cos(theta);
    
    exparg = xp.^2 + (aspectRatio^2)*yp.^2;
    cosarg = 2*pi*xp / freqBw;
    val = exp( -exparg / 2*tempBw^2 ) .* cos(cosarg + phase);
    