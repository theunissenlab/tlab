function [pxy, gridx, gridy] = kde2d_corrected(x, y, xbounds, ybounds)

    data = [rv(x) rv(y)];
    [bandwidth, pxy, X, Y, gridx, gridy] = kde2d(data, 2^7);
    
    minY = ybounds(1);
    maxY = ybounds(2);
    
    minX = xbounds(1);
    maxX = xbounds(2);    
        
    minIndxY = max(find(gridy < minY)) + 1;
    maxIndxY = min(find(gridy > maxY)) - 1;
    gridy = gridy(minIndxY:maxIndxY);
    
    minIndxX = max(find(gridx < minX)) + 1;
    maxIndxX = min(find(gridx > maxX)) - 1;
    gridx = gridx(minIndxX:maxIndxX);
    
    pxy = pxy(minIndxY:maxIndxY, minIndxX:maxIndxX);    
    pxy(pxy < 0) = 0;

    %% compute gaussian density map of original data
    maxX = max(data(:, 1));
    minX = min(data(:, 1));
    maxY = max(data(:, 2));
    minY = min(data(:, 2));

    xwidth = maxX - minX;
    ywidth = maxY - minY;
    
    spRadiusX = 0.3*xwidth;
    spRadiusY = 0.3*ywidth;
    
    [X, Y] = meshgrid(gridx, gridy);
    rspnts = [X(:) Y(:)];
    
    emap = zeros(size(pxy));
    
    for k = 1:size(rspnts, 1)        
        p = rspnts(k, :);
        xDens = density_in_rect(p, spRadiusX, 0.25*spRadiusY, data);
        yDens = density_in_rect(p, 0.25*spRadiusX, spRadiusY, data);
        
        if xDens > 0 && yDens > 0
            
            spx = 1/xDens;
            spy = 1/yDens;

            xcent = X - p(1);
            ycent = Y - p(2);

            g = (xcent.^2 / (2*spx^2)) + (ycent.^2 / (2*spy^2));
            emap = emap + exp(-g);        
        end
    end
    
    emap = emap / max(emap(:));
    
    %{
    figure; hold on;
    subplot(1, 4, 1); imagesc(pxy); axis tight; colorbar; title('before');
    subplot(1, 4, 2); imagesc(emap); axis tight; colorbar; title('idensity');
    
    pxyc = pxy .* emap;
    subplot(1, 4, 3); imagesc(pxyc); axis tight; colorbar; title('after');
    
    pxydiff = pxyc - pxy;
    subplot(1, 4, 4); imagesc(pxydiff); axis tight; colorbar; title('diff');
    %}
    pxyc = pxy;
    
    
end


function d = density_in_rect(center, width, height, data)

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
    
    d = n / width*height;    
end



    
    
    
