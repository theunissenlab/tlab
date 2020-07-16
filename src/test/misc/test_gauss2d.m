function test_gauss2d(bwx, bwy)

    x = -8:0.01:8;
    y = 0:0.05:1;
    
    cx = 0;
    cy = 0.5;
    
    [X, Y] = meshgrid(x, y);
    
    g = ((X-cx).^2 / (2*bwx)) + ((Y-cy).^2 / (2*bwy));
    g2d = exp(-g);
    
    figure; hold on;
    imagesc(x, y, g2d); axis tight; colorbar;
    title(sprintf('bx=%0.2f, by=%0.2f', bwx, bwy));
    
    