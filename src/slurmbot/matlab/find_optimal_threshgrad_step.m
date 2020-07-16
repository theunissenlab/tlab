function stepSize = find_optimal_threshgrad_step(modelParams, trainingIndex, validationIndex, optOptions, stepMult)

    if nargin < 5
        stepMult = 1e-4;
    end

    currentStep = 1;
    stepInc = 1;
    problemStep = 0;
    
    upTol = 1;
    
    while ~problemStep
        
        oo = optOptions;
        oo.maxIter = 20;
        oo.stepSize = currentStep*stepMult;
        
        fprintf('Testing out step size %f...\n', oo.stepSize);
        
        if optOptions.earlyStop
            [modelParamsOut, oo] = strfOpt(modelParams, trainingIndex, oo, validationIndex);
        else
            [modelParamsOut, oo] = strfOpt(modelParams, trainingIndex, oo);
        end
        
        errs = oo.diagnostics.errs(1, :);
        
        derrs = diff(errs);
        upCount = sum(derrs > 0);
        
        if upCount > upTol 
           break;
        end
        currentStep = currentStep + stepInc;
        
    end
    
    sdiff = 1;
    if currentStep > 2
        sdiff = 2;
    end
    
    stepSize = (currentStep - sdiff)*stepMult;
    fprintf('Optimal step size: %f\n', stepSize);
    
    