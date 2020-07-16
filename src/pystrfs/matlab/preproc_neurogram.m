function [ngram, groupIndex, ngramInfo] = preproc_neurogram(neurogramFile, stimMd5s, unitsToInclude, label)

    if nargin < 4
        label = 'best';
    end        

    h5 = h5utils();
    fid = h5.open(neurogramFile);
    
    nunits = length(unitsToInclude);
    
    stimLens = zeros(nunits, 1);
    allResponses = cell(nunits, 1);
    
    for k = 1:nunits
        resps = struct;        
        slen = 0;
        for j = 1:length(stimMd5s)            
            rpath = sprintf('/%s/%s/%s', stimMd5s{j}, unitsToInclude{k}, label);
            resp = h5.get_ds(fid, rpath);            
            md5Name = sprintf('md5_%s', stimMd5s{j});
            resps.(md5Name) = resp;
            %fprintf('(%s, %s), len=%d\n', unitsToInclude{k}, stimMd5s{j}, length(resp));
            stimLens(j) = length(resp); %should be same for all k
        end        
        allResponses{k} = resps;
    end
        
    h5.close(fid);
    
    %% construct neurogram
    totalLen = sum(stimLens);
    ngram = zeros(nunits, totalLen);    
    groupIndex = zeros(1, totalLen);
    s = 1;
    e = 0;
    for j = 1:length(stimMd5s)
 
        slen = stimLens(j);
        e = s + slen - 1;
        md5Name = sprintf('md5_%s', stimMd5s{j});
        for k = 1:nunits
            resp = allResponses{k}.(md5Name);
            ngram(k, s:e) = resp;
        end
        groupIndex(s:e) = j;
        s = e + 1;
    end
    
    %% normalize so all responses are between 0 and 1
    for k = 1:nunits
        ngram(k, :) = max(0, ngram(k, :));
        ngram(k, :) = ngram(k, :) / max(ngram(k, :));        
    end
    
    ngramInfo = struct;
    ngramInfo.stimLengths = stimLens;
    