function PlotSplit(moto, split)

    slope=-0.6025;
    intercept=-0.5979;

    mgram = moto.b;

    time = mgram(:, 1);    
    alpha = mgram(:, 2);  
    beta = slope * alpha + intercept;
    mu = mgram(:, 3);
    sigma1 = mgram(:, 4);
    sigma2 = mgram(:, 5);
    power = mgram(:, 6);
    fitvals = mgram(:, 7);

    findex = (split.f > 300) & (split.f < 10000);
    
    figure();    
    subplot(5, 1, 1); hold on;    
    imagesc(time, split.f(findex), split.spec(findex, :));
    set(gca,'YDir','normal');
    axis tight;
    ylabel('Frequency');
    
    subplot(5, 1, 2); hold on;    
    plot(time, alpha, 'b-');
    plot(time, beta, 'r-');
    axis tight;
    ylabel('alpha,beta');
    legend('alpha', 'beta');
    
    subplot(5, 1, 3); hold on;    
    plot(time, mu, 'm-');    
    ylabel('mu');
    axis tight;
    
    subplot(5, 1, 4); hold on;    
    plot(time, sigma1, 'r-');
    plot(time, sigma2, 'k-');
    legend('Sigma1', 'Sigma2');
    axis tight;
    ylabel('Sigma');    
    
    subplot(5, 1, 5); hold on;
    plot(time, fitvals, 'r-');
    axis tight;
    ylabel('Error');