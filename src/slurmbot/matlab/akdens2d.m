function pxy = akdens2d(data, xgrid, ygrid)

    maxX = max(data(:, 1));
    minX = min(data(:, 1));
    maxY = max(data(:, 2));
    minY = min(data(:, 2));

    xwidth = maxX - minX;
    ywidth = maxY - minY;
    
    npts = size(data, 1);

    spRadiusX = 0.05*xwidth;
    spRadiusY = 0.05*ywidth;
    
    [X, Y] = meshgrid(xgrid, ygrid);
    
    sparseness = zeros(npts, 2);
    
    %% compute sparseness for each data point
    for k = 1:npts       
        
        pnt = data(k, :);
        
        xCount = count_in_rect(pnt, spRadiusX, 0.25*spRadiusY, data);
        xDens = xCount / npts;
        
        yCount = count_in_rect(pnt, 0.25*spRadiusX, spRadiusY, data);
        yDens = yCount / npts;
        
        sparseness(k, :) = [xDens^-1 yDens^-1];
    end
    
    %% rescale sparseness so minimum encompasses minimum distance between
    %% data points and maximum encompasses maximum distance
    sdatax = sort(data(:, 1));
    dx = diff(sdatax);
    minPossSparseX = min(dx)
    %minPossSparseX = mean(dx)
    maxPossSparseX = max(dx)
    sxdist = maxPossSparseX - minPossSparseX;
    
    sdatay = sort(data(:, 2));
    dy = diff(sdatay);
    minPossSparseY = min(dy)
    %minPossSparseY = mean(dy)
    maxPossSparseY = max(dy)
    sydist = maxPossSparseY - minPossSparseY;
    
    sparsex = sparseness(:, 1);
    sparsey = sparseness(:, 2);
    
    sparsex = (sparsex/max(sparsex))*sxdist; %scale between min and max possible distances
    sparsey = (sparsey/max(sparsey))*sydist; %scale between min and max possible distances
    
    sparseness = [sparsex sparsey];
    
    sxmax = max(sparseness(:, 1))
    sxmin = min(sparseness(:, 1))
    
    symax = max(sparseness(:, 2))
    symin = min(sparseness(:, 2))
   
    %% compute probability map by summing gaussians across data points
    pxy = zeros(size(X));
    for k = 1:npts        
        xcent = X - data(k, 1);
        ycent = Y - data(k, 2);
        spx = sparseness(k, 1);
        spy = sparseness(k, 2);
        
        h = spx^2 + spy^2;
        g = (xcent.^2 / (2*spx^2)) + (ycent.^2 / (2*spy^2));
        pxy = pxy + h*exp(-g);
    end
    
end



function n = count_in_rect(center, width, height, data)

    hwidth = width / 2;
    hheight = height / 2;
    
    leftBound = center(1) - hwidth;
    rightBound = center(1) + hwidth;
    topBound = center(2) + hheight;
    bottomBound = center(2) - hheight;

    x = data(:, 1);
    y = data(:, 2);
    
    xindx = (x >= leftBound) & (x <= rightBound);
    yindx = (y >= bottomBound) & (y <= topBound);
    
    dindx = xindx & yindx;
    n = sum(dindx);
    
end

