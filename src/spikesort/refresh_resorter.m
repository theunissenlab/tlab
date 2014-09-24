function refresh_resorter(handles)
global partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels validtrials tmpv tmpv2 tmpclusters r5pc s5pc 
global heightsmean heightsstd widthsmean widthsstd goodflag fname tmpsortidx tmpunsortedidx clusize vs tridx tmps2 retry undohistory clu

bdf2 = disable_all(handles);
% 
% tmpsortidx(snpidx(find(labels==labelnum))) = clu;
% tmpunsortedidx = find(tmpsortidx==0);
% 
% if clu>0
% 
%     tempcluster
%     set(handles.figure1, 'Visible', 'on')
% end

% vv = cat(3,PF.rec(:).v);

kms = 6;

set(handles.ChanNum, 'String', ['Channel: ' num2str(ch)]);

undohistory{end+1}.tmpsortidx = tmpsortidx;
undohistory{end}.snpidx1 = snpidx1;
undohistory{end}.snpidx2 = snpidx2;
undohistory{end}.currsort=currsort;
undohistory{end}.partialflag = partialflag;
% undohistory{end+1}.tmpsortidx = tmpsortidx;

snpidx = intersect(snpidx1:snpidx2, tmpunsortedidx);
    


if ~isempty(tmpv2(:,snpidx))

    
    plot(handles.mainplot, tmpv2(:,snpidx)); set(handles.mainplot, 'YLim', [plotmin plotmax]);

    if get(handles.slope, 'Value')==1
        labels = cluster(diff(tmpv2(:,snpidx))', kms, handles);
    else
        labels = cluster(tmpv2(:,snpidx)', kms, handles);
    end


    for kk = 1:kms
        hh = ['handles.kmean' num2str(kk)];
        bdf = get(eval(hh), 'ButtonDownFcn');
        cla(eval(hh), 'reset')
        plot(eval(hh), tmpv2(:,snpidx(labels==kk))); set(eval(hh), 'YLim', [plotmin plotmax]);
        set(eval(hh), 'ButtonDownFcn', bdf)
    end
else


    if isempty(currsort)
        currsort=0;
    end

    flag = 0;
    kms = 6;

    while flag == 0

        if partialflag == 0 && ch < 97

            currsort = currsort + 1;

            if currsort > length(uvsorts)
                if goodflag==1
                    fc = finalcluster;
                    uiwait(fc)
                end
                
                goodflag = 0;
                
                if ~retry
                    ch = ch + 1;
                end
                
                if ch >= 97
                    continue
                end

                [vsorts, vtrials] = find(cellfun(@(x) ~isempty(x),vv(:,ch,:))==1);

                uvsorts = unique(vsorts);
                
                tmpv2 = zeros(size(cat(2,vv{uvsorts, ch, :})));
                tmps2 = zeros(1,size(tmpv2,2));
                tridx = [];
                uvidx = 0;
                clusize = [];
                
                uvsize=[];
        
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


                tmpsortidx = zeros(1, size(tmpv2,2));

                tmpunsortedidx = find(tmpsortidx==0);
                currsort = 1;
                if clusize(currsort) > str2num(get(handles.snipnum, 'string'));
                    snpidx1 = 1;
                    snpidx2 = str2num(get(handles.snipnum, 'string'));
                    partialflag = 1;
                else
                    snpidx1 = 1;
                    snpidx2 = clusize(currsort);
                    partialflag = 0;
                end
                
                undohistory = [];
                undohistory{1}.tmpsortidx=tmpsortidx;
                undohistory{1}.snpidx1=snpidx1;
                undohistory{1}.snpidx2=snpidx2;
                undohistory{1}.currsort=currsort;
                undohistory{1}.partialflag = partialflag;

            else

                snpidx1 = clusize(currsort-1)+1;
                if (clusize(currsort) - clusize(currsort-1)) > str2num(get(handles.snipnum, 'string'));
                    snpidx2 = snpidx1 + str2num(get(handles.snipnum, 'string'))-1;
                    partialflag = 1;
                else
                    if currsort == length(uvsorts)
                        snpidx2 = size(tmpv2,2);
                    else
                        snpidx2 = clusize(currsort);
                    end
                    partialflag = 0;
                end
                
                if isempty(intersect(snpidx1:snpidx2, tmpunsortedidx))
                    continue
                end

            end

            set(handles.ChanNum, 'String', ['Channel: ' num2str(ch)]);



        elseif partialflag == 1 && ch < 97

            snpidx1 = snpidx1 + str2num(get(handles.snipnum, 'string'));
            snpidx2 = snpidx2 + str2num(get(handles.snipnum, 'string'));

            if clusize(currsort) < snpidx2
                snpidx2 = clusize(currsort);
                partialflag = 0;
            else
                partialflag = 1;
            end
        
        else
%             load(fname)
             s1 = strfind(fname, '/');
% 
%             for vt = validtrials
%                 PF.rec(vt).v_final = tmpclusters(vt).v;
%                 PF.rec(vt).plexon_spikes_final = tmpclusters(vt).plexon_spikes;
%             end

            save([fname(1:end-3) 'rs.mat'], 'PF')

            save([fname(1:s1(end)) 'sortparams.mat'], 'ens', 'heightsmean', 'heightsstd', 'widthsmean', 'widthsstd', 'r5pc', 's5pc')
            fprintf('Done')
            return
            
        end


%         snpidx = snpidx1:snpidx2;
         snpidx = intersect(snpidx1:snpidx2, tmpunsortedidx);

%         tmpv2 = cat(2,vv{uvsorts, ch, :});

%         if mean(tmpv2(1,:))<mean(tmpv2(7,:))
% 
%             continue
%         else
            goodflag=1;
%         end

        if (currsort == 1 && snpidx1 == 1) || isempty(heightsmean)


            [plotmin, plotmax] = findplotminmax(tmpv2);
            allclusters
%             tmpv2 = cat(2,vv{uvsorts, ch, :});
            
            covtv = diff(tmpv2)*diff(tmpv2)';
            [U,S,V]= svd(covtv);

            s5pc{ch} = U;

            covtv = tmpv2*tmpv2';
            [U,S,V]= svd(covtv);

            r5pc{ch} = U;

            heights = max(tmpv2) - min(tmpv2);
            heightsmean = mean(heights);
            heights = heights - heightsmean;
            heightsstd = std(heights);
            [minr,minc] = find(bsxfun(@(x,y) (x == y), tmpv2, min(tmpv2))==1);
            [maxr,maxc] = find(bsxfun(@(x,y) (x == y), tmpv2, max(tmpv2))==1);
            minr = minr(unique(minc));
            maxr = maxr(unique(maxc));
            widths = minr - maxr;
            widthsmean = mean(widths);
            widths = widths - widthsmean;
            widthsstd = std(widths);

            tempcluster
            set(handles.figure1, 'Visible', 'on')
        end


        plot(handles.mainplot, tmpv2(:,snpidx)); set(handles.mainplot, 'YLim', [plotmin plotmax]);

        if get(handles.slope, 'Value')==1
            labels = cluster(diff(tmpv2(:,snpidx))', kms, handles);
        else
            labels = cluster(tmpv2(:,snpidx)', kms, handles);
        end

        for kk = 1:kms
            hh = ['handles.kmean' num2str(kk)];
            bdf = get(eval(hh), 'ButtonDownFcn');
            cla(eval(hh), 'reset')
            plot(eval(hh), tmpv2(:,snpidx(labels==kk))); set(eval(hh), 'YLim', [plotmin plotmax]);
            set(eval(hh), 'ButtonDownFcn', bdf)
        end

        flag = 1;
    end
end
enable_all(handles, bdf2)

function bdf = disable_all(handles)
set(handles.status, 'String', 'Wait');  drawnow
set(handles.nextbutton, 'Enable', 'off')
set(handles.skip, 'Enable', 'off')
bdf{1} = get(handles.kmean1, 'ButtonDownFcn');
bdf{2} = get(handles.kmean2, 'ButtonDownFcn');
bdf{3} = get(handles.kmean3, 'ButtonDownFcn');
bdf{4} = get(handles.kmean4, 'ButtonDownFcn');
bdf{5} = get(handles.kmean5, 'ButtonDownFcn');
bdf{6} = get(handles.kmean6, 'ButtonDownFcn');
set(handles.kmean1, 'ButtonDownFcn', []);
set(handles.kmean2, 'ButtonDownFcn', []);
set(handles.kmean3, 'ButtonDownFcn', []);
set(handles.kmean4, 'ButtonDownFcn', []);
set(handles.kmean5, 'ButtonDownFcn', []);
set(handles.kmean6, 'ButtonDownFcn', []);

function enable_all(handles, bdf)
set(handles.status, 'String', 'Ready');
set(handles.nextbutton, 'Enable', 'on');
set(handles.skip, 'Enable', 'on');
set(handles.kmean1, 'ButtonDownFcn', bdf{1});
set(handles.kmean2, 'ButtonDownFcn', bdf{2});
set(handles.kmean3, 'ButtonDownFcn', bdf{3});
set(handles.kmean4, 'ButtonDownFcn', bdf{4});
set(handles.kmean5, 'ButtonDownFcn', bdf{5});
set(handles.kmean6, 'ButtonDownFcn', bdf{6});

