function test_glogistic()
    
    clrs = {'k-', 'b-', 'c-', 'g-', 'y-', 'r-'};

    Q = [0.001 0.01 0.1 1 5];
    B = [0.001 0.01 0.1 1 2 5];
    M = [0.001 0.01 0.1 1 2 5];
    
    x = -10:0.05:10;
    
    figure; hold on;
    for k = 1:length(Q)
        qval = Q(k);
        vval = Q(k);
        bval = 1;
        mval = 1;
        y = glogistic(x, qval, bval, mval, vval);
        plot(x, y, clrs{k});
    end
    title('Q');
    legend('0.001', '0.01', '0.1', '1', '5');
    axis tight;
    
    
    figure; hold on;
    for k = 1:length(B)
        bval = B(k);
        vval = 1;
        qval = 1;
        mval = 1;
        y = glogistic(x, qval, bval, mval, vval);
        plot(x, y, clrs{k});
    end
    title('B');
    legend('.001', '.01', '.1', '1', '2', '5');
    axis tight;
    
    figure; hold on;
    for k = 1:length(M)
        bval = 1;
        vval = 1;
        qval = 1;
        mval = M(k);
        y = glogistic(x, qval, bval, mval, vval);
        plot(x, y, clrs{k});
    end
    title('M');
    legend('.001', '.01', '.1', '1', '2', '5');
    axis tight;
    
    