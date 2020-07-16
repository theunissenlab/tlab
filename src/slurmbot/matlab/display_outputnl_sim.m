function display_outputnl_sim(ofile)

    %vars = load(ofile);
    load(ofile);    
        
    %% plot simulation
    figure; hold on;
    subplot(2, 1, 1); hold on;
    plot(fullResponse.stimCurrent, 'g-');
    axis tight;
    subplot(2, 1, 2); hold on;
    plot(modelResponse, 'k-');
    axis tight;
    
    oftypes = {'err', 'kl'};
    
    for k = 1:length(oftypes)        
        oftype = oftypes{k};
        psthTrain = perfInfo{1}.psthTrain;
        psthValid = perfInfo{1}.psthValid;
        
        perfInfoDists = perfInfo{1}.(oftype);
        perfInfoSpline = perfInfo{2}.(oftype);
        nlinfoDists = perfInfo{1}.nlinfo.(oftype);
        nlinfoSpline = perfInfo{2}.nlinfo.(oftype);
    
        %% plot distributions and output nonlinearities    
        figure('Name', oftype); hold on;

        subplot(2, 1, 1); hold on;
        plot(nlinfoDists.x, nlinfoDists.px_smooth, 'k-'); axis tight;
        plot(nlinfoDists.x, nlinfoDists.pxspike_smooth, 'r-'); axis tight;
        plot(nlinfoDists.x, nlinfoDists.pxnospike_smooth, 'b-'); axis tight;    
        legend('p(x)', 'p(x|spike)', 'p(x|nospike)');

        subplot(2, 1, 2); hold on;

        splineOut = fnval(nlinfoSpline.outputNL, nlinfoDists.x);
        distOut = fnval(nlinfoDists.outputNL, nlinfoDists.x);

        axval = [min(nlinfoDists.x) max(nlinfoDists.x) 0 2];

        plot(nlinfoDists.x, nlinfoDists.pspikex_smooth, 'r-');
        plot(nlinfoDists.x, nlinfoDists.pnospikex_smooth, 'b-');
        plot(nlinfoDists.x, distOut, 'm-');
        plot(nlinfoDists.x, splineOut, 'g-');
        axis(axval);
        legend('p(spike|x)', 'p(nospike|x)', '~p(spike|x)', 'spline');
        
        
        %% plot nonlinear responses
        figure('Name', oftype); hold on;    

        subplot(2, 1, 1); hold on;
        plot(psthValid, 'k-');
        plot(perfInfoDists.validResp, 'b-');
        axis tight;
        title(sprintf('NL from p(spike|x): perf=%0.3f', perfInfoDists.perfRatioValid));

        subplot(2, 1, 2); hold on;
        plot(psthValid, 'k-');
        plot(perfInfoSpline.validResp, 'b-');
        axis tight;
        title(sprintf('NL from spline: perf=%0.3f', perfInfoSpline.perfRatioValid));

        fprintf('%s: NL from p(spike|x): train=%0.5f  |  valid=%0.5f\n', oftype, perfInfoDists.perfRatioTrain, perfInfoDists.perfRatioValid);
        fprintf('%s: NL from spline: train=%0.5f  |  valid=%0.5f\n', oftype, perfInfoSpline.perfRatioTrain, perfInfoSpline.perfRatioValid);

    end
    