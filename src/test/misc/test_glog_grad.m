function test_glog_grad()

    M = 0.5;
    B = 1.2;
    
    xvals = -5:0.1:5;
    yvals = glogistic(xvals, B, M);
    %yvals = randn(size(yvals))*1e-3 + yvals;
    
    datIdx = 1:length(yvals);
    groupIndex = ones(size(yvals));
    
    strfData(xvals, yvals, groupIndex);
    
    gp = glogInit();
    gp.M = M;
    gp.B = B;
    
    %test glogFwd    
    %{
    [gp, presp] = glogFwd(gp, datIdx);    
    figure; hold on;
    plot(xvals, yvals, 'k-');
    plot(xvals, presp, 'r-');
    %}
    
    
    gp.M = randn();
    gp.B = abs(randn())*1.25;
    
    [gp, g] = glogGrad(gp, datIdx);    
    [gp, gfd] = glogGradFD(gp, datIdx);
    
    g
    gfd
    norm(g - gfd)
    