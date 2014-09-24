function [avg_strf, avg_bias, avg_output_nl] = compute_avg_strf_and_nl(responseFileTemplate, nfolds)

    if nargin < 2
        nfolds = 5;
    end
    
    strf_list = cell(nfolds, 1);
    biases = zeros(nfolds, 1);
    outputNLs = cell(nfolds, 1);
    
    numChannels = -1;
    numDelays = -1;
       
    max_minx = -Inf;
    min_maxx = Inf;
    
    for k = 1:nfolds
    
        responseFile = sprintf(responseFileTemplate, k-1);
        [strf, bias, outputNL, modelType, x, m] = get_strf_and_nl(responseFile);
        
        minx = min(x);
        maxx = max(x);
        
        if minx > max_minx
            max_minx = minx;
        end
        if maxx < min_maxx
            min_maxx = maxx;
        end
        
        sz = size(strf);
        if numChannels == -1
            numChannels = sz(1);
        end
        if numDelays == -1
            numDelays = sz(2);
        end
       
        strf_list{k} = strf;
        biases(k) = bias;
        if ~strcmp(modelType, 'leglm')
            outputNLs{k} = outputNL;
        else
            outputNLs{k} = m;
        end
    end
    
    strfs = zeros(numChannels, numDelays, nfolds);
    for k = 1:nfolds       
        strfs(:, :, k) = strf_list{k};
    end
    clear strf_list;
    
    avg_strf = mean(strfs, 3);
    avg_bias = mean(biases);
    if ~ismember(modelType, {'nl_dists', 'nl_spline', 'leglm'})
        avg_output_nl = @(x) outputNLs{1}(x);
    else
        if strcmp(modelType, 'leglm')
            avg_m = mean(cell2mat(outputNLs))
            avg_output_nl = @(x) log(1 + exp(x)).^avg_m; 
        else
            avg_output_nl = @(x) mean_output_nl(x, outputNLs, max_minx, min_maxx);
        end
    end
    
end


function [strf, bias, outputNL,modelType, x, m] = get_strf_and_nl(responseFile)

    h5 = h5utils();
    fid = h5.open(responseFile);
    
    m = -1;

    modelType = h5.get_attr(fid, '/model', 'type');
       
    %% get STRF    
    if ~ismember(modelType, {'nl_dists', 'nl_spline'})        
        strf = h5.get_ds(fid, '/model/weights');
        bias = h5.get_ds(fid, '/model/bias');
    else
        origFile = h5.get_attr(fid, '/model', 'original_response_file');
        fid2 = h5.open(origFile);
        strf = h5.get_ds(fid2, '/model/weights');
        bias = h5.get_ds(fid2, '/model/bias');
        h5.close(fid2);
    end
    
    x = h5.get_ds(fid, '/model/output_nl/domain');
    
    %% get output NL
    if ismember(modelType, {'nl_dists', 'nl_spline'})
        [rootDir, baseName, ext] = fileparts(responseFile);
        matFile = fullfile(rootDir, [baseName '.mat']);
        vars = load(matFile);
        outputNL = @(x) output_nl_wrapper(x, vars.nlinfo.outputNL);
        clear vars;
    else       
        switch modelType            
            case 'poisson'
                outputNL = @(x) exp(x);                
            case 'binomial'
                outputNL = @(x) 1 ./ (1 + exp(-x));                
            case 'leglm'
                m = h5.get_ds(fid, '/model/m');
                outputNL = @(x) log(1 + exp(x)).^m;                
            case 'linear'
                outputNL = @(x) x;
        end
    end
    
    h5.close(fid);

end

function y = mean_output_nl(x, meanOutputNLs, max_minx, min_maxx)

    n = length(meanOutputNLs);
    y = zeros(size(x));
    for k = 1:n
        y = y + meanOutputNLs{k}(x);
    end
    y(x < max_minx) = 0;
    y(x > min_maxx) = 0;
    y = y / n;
end

function [y, dy] = output_nl_wrapper(x, fnObj)    
    y = fnval(fnObj, x);
    if nargout > 1
        dy = zeros(size(y));
    end
end

