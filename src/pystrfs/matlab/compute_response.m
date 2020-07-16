function [resp, rawResp] = compute_response(stimRespData, modelParamsTrained, rectify)

    if nargin < 3
       rectify = 1;
    end

    wholeResp = stimRespData.wholeResp;
    
    [modelParamsTrained, rawResp] = strfFwd(modelParamsTrained, 1:length(wholeResp));
    
    rawResp(isnan(rawResp)) = 0;
    resp = rawResp;
    if rectify
        resp(resp < 0) = 0;
        resp = (resp / max(resp))*max(wholeResp);
    end
    