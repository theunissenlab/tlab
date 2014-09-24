function fig_glm_outputnls()


    %{
    figure; hold on;
    x = -3:0.1:3;
    plot(x, exp(x), 'k-', 'LineWidth', 6);
    axis tight;
    
    figure; hold on;
    x = -5:0.1:5;
    plot(x, 1 ./ (1 + exp(-x)), 'k-', 'LineWidth', 6);
    axis tight;
    %}

    %{
    prand = poissrnd(10, 1000, 1);
    figure; hold on;
    hist(prand, 10);
    %}

    brand = binornd(20, 0.5, 5000, 1);
    figure; hold on;
    hist(brand, 20);
    

    %{
    nrand = randn(50000, 1);
    figure; hold on;
    hist(nrand, 75);
    %}


    %{
    figure; hold on;
    x = -5:0.1:5;
    mvals = [1 1.5 2];
    clrs = {'b-', 'g-', 'r-'};
    for k = 1:length(mvals)
        m = mvals(k);
        plot(x, log(1 + exp(x)).^m, clrs{k}, 'LineWidth', 6);
    end
    axis tight;
    %}

    %{
    figure; hold on;
    x = -5:0.1:5;
    plot(x, x, 'k-', 'LineWidth', 6);
    %}
