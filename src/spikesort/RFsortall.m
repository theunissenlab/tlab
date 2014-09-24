function RFsortall(filenum)

javaaddpath('/usr/share/java/mysql-connector-java.jar');
addpath('/auto/k1/shinji/matlab/arraydb')
d=arraydbgettask('runclass','FP30');


fpname=[d(filenum).rawtaskpath d(filenum).rawtaskfile '.tdt.rs.mat']

if exist(fpname) ~= 2
    error('Not yet sorted')
    return
end

fp = load(fpname);

vvfp = cat(3,fp.PF.rec(:).v_final);

for ii = 1:size(fp.PF.rec,2)
    vsizes(ii) = size(fp.PF.rec(ii).v_final,1);
end

vvv = find(vsizes == max(vsizes));
vvv= vvv(1);


sortfilename = [d(filenum).rawtaskpath 'sortparams.mat'];

sortparams = load(sortfilename);

sortedchans = find(cellfun(@(x) ~isempty(x), sortparams.ens)==1);

a = dir([d(filenum).rawtaskpath '*tdt.mat']);

for aa = 1:size(a,1)

    sortfiles{aa} = [d(filenum).rawtaskpath a(aa).name];
end

for ss = 1:size(sortfiles,2)
    
    outfname = [sortfiles{ss}(1:end-3) 'rs.mat'];
    
    if exist(outfname)~=2

        PFname = sortfiles{ss};

        clear PF vsizes
        
        load(PFname)

        for ii = 1:size(PF.rec,2)
            if ~isfield(PF.rec(ii), 'v')
                vsizes(ii)=0;
            else
                vsizes(ii) = size(PF.rec(ii).v,1);
            end
        end
        
        if max(vsizes) > 0 

        fixidx = find(vsizes < max(vsizes));

        for ii = fixidx
            if size(PF.rec(ii).v,1)>0
                PF.rec(ii).plexon_spikes{max(vsizes),1}=[];
                PF.rec(ii).v{max(vsizes),1} = [];
            end
        end

        for ii = 1:size(PF.rec,2)
            vsizes(ii) = size(PF.rec(ii).v,2);
        end

        fixidx = find(vsizes < max(vsizes));

        for ii = fixidx
            if size(PF.rec(ii).v,1)>0
                PF.rec(ii).plexon_spikes{1,max(vsizes)}=[];
                PF.rec(ii).v{1,max(vsizes)} = [];
            end
        end

        vv = cat(3,PF.rec(:).v);
        vs = cat(3,PF.rec(:).plexon_spikes);
        vt = 1:size(PF.rec,2);
        validtrials=[];
        for ii = vt
            if isempty(PF.rec(ii).v)
                good = 0;
            else
                good = 1;
            end
            validtrials = [validtrials good];
        end

        validtrials = vt(logical(validtrials));

        for vt = validtrials
                    PF.rec(vt).v_final = cell(size(fp.PF.rec(vvv).v_final));
                    PF.rec(vt).plexon_spikes_final = cell(size(fp.PF.rec(vvv).plexon_spikes_final));
                    PF.rec(vt).F_final = cell(size(fp.PF.rec(vvv).F_final));
                    PF.rec(vt).type_final = cell(size(fp.PF.rec(vvv).type_final));
        end

        for ch = sortedchans

            [vsorts, vtrials] = find(cellfun(@(x) ~isempty(x),vv(:,ch,:))==1);

            uvsorts = unique(vsorts);

            tmpv2 = zeros(size(cat(2,vv{uvsorts, ch, :})));
            tmps2 = zeros(1,size(tmpv2,2));
            tridx = [];
            uvidx = 0;
            clusize = [];
            uvsize=[];
            
            if ~isempty(tmpv2)

                for uv = uvsorts'
                    uvsize = [uvsize size(cat(2,vv{uv, ch, :}),2)];
                end

                [~,uvs] = sort(uvsize, 'descend');
                uvsorts = uvsorts(uvs);

                for uv = uvsorts'
                    tmpv = cat(2,vv{uv, ch, :});
                    tmpv2(:,uvidx+1:uvidx+size(tmpv,2)) = tmpv;
                    tmps = cat(2,vs{uv, ch, :});
                    tmps2(:,uvidx+1:uvidx+size(tmps,2)) = tmps;
                    uvidx = uvidx+size(tmpv,2);
                    clusize = [clusize uvidx];

                    spintr = cat(1,squeeze(cellfun(@(x) size(x,2), vv(uv, ch, :))));

                    for ii = 1:length(spintr)
                          tridx = [tridx; repmat(validtrials(ii), spintr(ii),1)];
                    end
                end
                %         keyboard

                tmpsortidx = zeros(1, size(tmpv2,2));
                tmpunsortedidx = find(tmpsortidx==0);


                if strcmp(class(sortparams.ens{ch}), 'double') && sortparams.ens{ch} == 1

                    labels = ones(size(tmpv2(:,tmpunsortedidx),2),1);
                    kms=1;
                else
                    X = zeros(101, size(tmpv2(:,tmpunsortedidx),2));
                    X(1:30,:) = tmpv2(:,tmpunsortedidx);
                    X(31:59,:) = diff(tmpv2(:,tmpunsortedidx));
                    X(60:89,:) = abs(fft(tmpv2(:,tmpunsortedidx)));

                    if ch > length(sortparams.s5pc)

                        tmpvfp = cat(2,vvfp{:, ch, :});

                        covtv = diff(tmpvfp)*diff(tmpvfp)';
                        [U,S,V]= svd(covtv);

                        sortparams.s5pc{ch} = U;

                        covtv = tmpvfp*tmpvfp';
                        [U,S,V]= svd(covtv);

                        sortparams.r5pc{ch} = U;

                    end

                    X(90:94,:) = sortparams.r5pc{ch}(:,1:5)'*tmpv2(:,tmpunsortedidx);
                    X(95:99,:) = sortparams.s5pc{ch}(:,1:5)'*X(31:59,:);
                    heights = max(tmpv2) - min(tmpv2(:,tmpunsortedidx));
                    X(100,:) = heights;
                    [minr,minc] = find(bsxfun(@(x,y) (x == y), tmpv2(:,tmpunsortedidx), min(tmpv2(:,tmpunsortedidx)))==1);
                    [maxr,maxc] = find(bsxfun(@(x,y) (x == y), tmpv2(:,tmpunsortedidx), max(tmpv2(:,tmpunsortedidx)))==1);
                    minr = minr(unique(minc));
                    maxr = maxr(unique(maxc));
                    widths = minr - maxr;
                    X(101,:) = widths;

                    labels = predict(sortparams.ens{ch}, X');
                    labels = str2num(cell2mat(labels));
                end

                tmpsortidx= labels';

                for vt = validtrials
                    for ll = unique(labels)'
                        PF.rec(vt).v_final{ll,ch} = tmpv2(:,((tridx== vt) &  (tmpsortidx' == ll)));
                        PF.rec(vt).plexon_spikes_final{ll,ch} = tmps2(((tridx== vt) &  (tmpsortidx' == ll)));
                        if ~isempty(tmpv2(:,((tridx== vt) &  (tmpsortidx' == ll))))
                            if ch > size(fp.PF.rec(vvv).F_final,2) || ll > size(fp.PF.rec(vvv).F_final,1)
                                PF.rec(vt).F_final{ll,ch} = [];
                            else
                                PF.rec(vt).F_final{ll,ch} = fp.PF.rec(vvv).F_final{ll,ch};
                            end
                            if ch > size(fp.PF.rec(vvv).F_final,2) || ll > size(fp.PF.rec(vvv).F_final,1)
                                PF.rec(vt).F_final{ll,ch} = 'multi';
                            elseif strcmp(fp.PF.rec(vvv).type_final{ll,ch}, 'single');
                                PF.rec(vt).type_final{ll,ch} = 'single';
                            elseif strcmp(fp.PF.rec(vvv).type_final{ll,ch},'multi');
                                PF.rec(vt).type_final{ll,ch} = 'multi';
                            elseif strcmp(fp.PF.rec(vvv).type_final{ll,ch}, 'noise');
                                PF.rec(vt).type_final{ll,ch} = 'noise';
                            end
                        end
                    end
                end
            end
            end
        
        save(outfname, 'PF')
        fprintf(['Saved ' outfname '\n'])
        end
    end
end


