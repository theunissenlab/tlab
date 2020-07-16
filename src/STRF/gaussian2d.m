function val = gaussian2d(x, y, x0, y0, params)
   
    xspread = params(1);
    yspread = params(2);
    xpart = (x - x0).^2 / 2*xspread^2;
    ypart = (y - y0).^2 / 2*yspread^2;
    val = exp(-xpart - ypart);
    
        