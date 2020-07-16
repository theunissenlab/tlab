function test_kernel2d()

    data = randn(2000, 2);
        
    [bandwidth,density,X,Y,xvals,yvals] = kde2d(data, 2^6);    
    
    xgrid = -4:0.05:4;
    ygrid = xgrid;
    %pxy = akdens2d(data, xgrid, ygrid);
    pxy = kde2d_nn(data, xgrid, ygrid);
    
    figure; hold on;
    plot(data(:, 1), data(:, 2), 'kx');
    axis([-4 4 -4 4]);
    
    figure; hold on;
    imagesc(xgrid, ygrid, pxy); axis tight; colorbar;
    title('akdens2d');
    
    figure; hold on;
    imagesc(xvals, yvals, density); axis tight; colorbar;
    title('kde2d');
    
    