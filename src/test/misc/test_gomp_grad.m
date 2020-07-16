function test_gomp_grad()

    a = 1;
    b = -5;
    c = -2;

    xvals = -5:0.1:5;
    yvals = gompertz(a, b, c, xvals);
    yvals = randn(size(yvals))*1e-2 + yvals;
    
    datIdx = 1:length(yvals);
    groupIndex = ones(size(yvals));
    
    strfData(xvals, yvals, groupIndex);
    
    gp = gompInit(a);
    %gp.b = b;
    %gp.c = c;
    
    %{
    %test gompFwd
    gp.b = b;
    gp.c = c;
    [gp, presp] = gompFwd(gp, datIdx);    
    figure; hold on;
    plot(xvals, yvals, 'k-');
    plot(xvals, presp, 'r-');
    %}
    
    %xpnt = randn(1, 2);
    %gp.b = -abs(xpnt(1)*2);
    %gp.c = -abs(xpnt(2)*2);
    gp.b = -2;
    gp.c = -0.5;
    
    [gp, g] = gompGrad(gp, datIdx);    
    [gp, gfd] = gompGradFD(gp, datIdx);
    
    g
    gfd
    norm(g - gfd)
    
    