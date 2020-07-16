function test_spikedists_cv()

    %mdata = get_model_data('/auto/k6/mschachter/pystrfs/units/pupi0414_10/output/directfit.Con.stft.nstd_6.fband_125.h5');
    mdata = get_model_data('/auto/k6/mschachter/pystrfs/units/yy1617_5_B/output/directfit.Con.stft.nstd_6.fband_125.h5');
    
    linResp = mdata.model.rawResponse;
    psth = mdata.data.wholeResp;    
    numTrials = round(mean(mdata.data.numTrials));
    
    numDensityPoints = 150;
    numResample = 5000;
    resampleFraction = 0.75;
    
    distinfo = estimate_spike_dists_cv(linResp, psth, numTrials, numResample, resampleFraction);
    %distinfo = estimate_spike_dists_cv2(linResp, psth, numTrials, numDensityPoints, numResample, resampleFraction);

    x = distinfo.x;
    
    figure; hold on;    
    subplot(3, 1, 1); hold on;
    hist(linResp, numDensityPoints);
    
    subplot(3, 1, 2); hold on;
    plot(x, distinfo.px, 'k-');
    plot(x, distinfo.pxspike, 'r-');
    plot(x, distinfo.pxnospike, 'b-');
    legend('p(x)', 'p(x|spike)', 'p(x|nospike)');
    axis tight;
    
    subplot(3, 1, 3); hold on;
    plot(x, distinfo.pspikex, 'r-', 'LineWidth', 2);
    plot(x, distinfo.pnospikex, 'b-');
    legend('p(spike|x)', 'p(nospike|x)');
    axis tight;    
    
    figure; hold on;
    plot(x, fnval(distinfo.outputNL, x), 'k-');
    plot(x, distinfo.pspikex, 'r-');
    title('Best Fit');
    axis tight;    
    