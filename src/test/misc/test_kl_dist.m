function test_kl_dist()

    npnts = 100000;
    
    a = 0;
    b = 0.0;
    
    x = randn(npnts, 1) + a;
    y = randn(npnts, 1) + b;
    
    dx = 0.001;
    rng = -15:dx:15;
    
    [px, xi, bwpx] = ksdensity(x, rng);
    [py, yi, bwpy] = ksdensity(y, rng);
    
    d = kl_dist(px, py, rng)    
    
    
    p1 = exp(-x.^2 / 2);
    c = (sqrt(2*pi)*2*log2(2));
    dpred = (sum(p1 .* (b^2 - 2*x*b)) / c)*dx
    
    figure; hold on;
    plot(px, 'r-');
    plot(py, 'b-');
    legend('p(x)', 'p(y)');
    
    