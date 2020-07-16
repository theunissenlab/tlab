function varargout = finalcluster_single(varargin)
% FINALCLUSTER_SINGLE MATLAB code for finalcluster_single.fig
%      FINALCLUSTER_SINGLE, by itself, creates a new FINALCLUSTER_SINGLE or raises the existing
%      singleton*.
%
%      H = FINALCLUSTER_SINGLE returns the handle to a new FINALCLUSTER_SINGLE or the handle to
%      the existing singleton*.
%
%      FINALCLUSTER_SINGLE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FINALCLUSTER_SINGLE.M with the given input arguments.
%
%      FINALCLUSTER_SINGLE('Property','Value',...) creates a new FINALCLUSTER_SINGLE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before finalcluster_single_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to finalcluster_single_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help finalcluster_single

% Last Modified by GUIDE v2.5 01-May-2012 14:29:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @finalcluster_single_OpeningFcn, ...
                   'gui_OutputFcn',  @finalcluster_single_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before finalcluster_single is made visible.
function finalcluster_single_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to finalcluster_single (see VARARGIN)
global PF vv ch currsort uvsorts vtrials vsorts labels sort validtrials tmpv2 tmpclusters sortparams kms kmslast plotmin plotmax tmpunsortedidx tmpsortidx vs tridx tmps2
% Choose default command line output for finalcluster_single
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% 
% UIWAIT makes finalcluster_single wait for user response (see UIRESUME)
% uiwait(handles.figure1);

if ~isempty(sortparams.ens{ch})
    
        sortedchans = find(cellfun(@(x) ~isempty(x), sortparams.ens)==1);
        
        scidx = find(sortedchans ==ch);
        
        if scidx+1 > length(sortedchans)
            return
        end
        
        if isempty(scidx)
            ch = sortedchans(1);
        else
            ch = sortedchans(scidx+1);
        end
        
        [vsorts, vtrials] = find(cellfun(@(x) ~isempty(x),vv(:,ch,:))==1);

        uvsorts = unique(vsorts);

        tmpv2 = zeros(size(cat(2,vv{uvsorts, ch, :})));
        tmps2 = zeros(1,size(tmpv2,2));
        tridx = [];
        uvidx = 0;
        clusize = [];
        
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
%         

        tmpsortidx = zeros(1, size(tmpv2,2));
        tmpunsortedidx = find(tmpsortidx==0);
        
        [plotmin, plotmax] = findplotminmax(tmpv2)

        
        if strcmp(class(sortparams.ens{ch}), 'double') && sortparams.ens{ch} == 1
            
            labels = ones(size(tmpv2(:,tmpunsortedidx),2),1);
            kms=1;
        else
            X = zeros(101, size(tmpv2(:,tmpunsortedidx),2));
            X(1:30,:) = tmpv2(:,tmpunsortedidx);
            X(31:59,:) = diff(tmpv2(:,tmpunsortedidx));
            X(60:89,:) = abs(fft(tmpv2(:,tmpunsortedidx)));
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

            kms = str2num(sortparams.ens{ch}.ClassNames{end});
        end
        
        set(handles.numcluster, 'String', num2str(kms));
else


    kms = max(1, length(unique(tmpsortidx(tmpsortidx>0))));

    set(handles.numcluster, 'String', num2str(kms));

    labels = cluster(tmpv2(:,tmpunsortedidx), kms, handles);

end

for ll = unique(labels)'
    tmpsortidx(tmpunsortedidx(labels==ll)) = ll;
end

kmslast = kms;
set(handles.ChanNum, 'String', ['Channel: ' num2str(ch)]);

    for kk = 1:6
        hh = ['handles.axes' num2str(kk)];
        bdf = get(eval(hh), 'ButtonDownFcn');
        cla(eval(hh), 'reset')
        plot(eval(hh), tmpv2(:,tmpunsortedidx(labels==kk))); set(eval(hh), 'YLim', [plotmin plotmax]);
        hold(eval(hh)); plot(eval(hh), mean(tmpv2(:,tmpunsortedidx(labels==kk)),2), 'k', 'linewidth', 4)
        set(eval(hh), 'ButtonDownFcn', bdf)
    end

    if strcmp(PF.extradata{2}.data, 'fp30.index')
        scorefunc(handles)
    end
        
% --- Outputs from this function are returned to the command line.
function varargout = finalcluster_single_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function numcluster_Callback(hObject, eventdata, handles)
% hObject    handle to numcluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numcluster as text
%        str2double(get(hObject,'String')) returns contents of numcluster as a double


% --- Executes during object creation, after setting all properties.
function numcluster_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numcluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function numpcs_Callback(hObject, eventdata, handles)
% hObject    handle to numpcs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numpcs as text
%        str2double(get(hObject,'String')) returns contents of numpcs as a double


% --- Executes during object creation, after setting all properties.
function numpcs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numpcs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in kmeans.
function kmeans_Callback(hObject, eventdata, handles)
% hObject    handle to kmeans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global PF vv ch currsort uvsorts vtrials vsorts labels sort validtrials tmpv tmpv2 subspace plotmin plotmax kmslast tmpunsortedidx tmpsortidx
set(handles.status, 'String', 'Wait'); drawnow
kms = str2num(get(handles.numcluster, 'String'));

if get(handles.randomforrest, 'Value') == 1
    set(handles.kmeanscluster, 'Value', 1)
    set(handles.spectralcluster, 'Value',0)
    set(handles.kmedoidscluster, 'Value',0)
    set(handles.randomforrest, 'Value',0)
end

subspace = 0;
% tmpv = cat(2,vv{uvsorts, ch, :});

if kms ~= kmslast
    tmpsortidx = zeros(1, size(tmpv2,2));
    tmpunsortedidx = find(tmpsortidx==0);
    kmslast = kms;
end

if kms == 1
    labels = ones(size(tmpv2,2),1);
else

    labels = cluster(tmpv2(:,tmpunsortedidx)', kms, handles);
end

for ll = unique(labels)'
    tmpsortidx(tmpunsortedidx(labels==ll)) = ll;
end

for kk = 1:6
    hh = ['handles.axes' num2str(kk)];
    bdf = get(eval(hh), 'ButtonDownFcn');
    cla(eval(hh), 'reset')
    plot(eval(hh), tmpv2(:,(labels==kk))); hold(eval(hh)); plot(eval(hh),mean(tmpv2(:,(labels==kk)),2), 'k', 'linewidth', 4);
    set(eval(hh), 'YLim', [plotmin plotmax]);
    set(eval(hh), 'ButtonDownFcn', bdf)
end
if strcmp(PF.extradata{2}.data, 'fp30.index')
    scorefunc(handles)
end
 set(handles.status, 'String', 'Ready'); drawnow
% --- Executes on button press in pca.
function pca_Callback(hObject, eventdata, handles)
% hObject    handle to pca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global PF vv ch currsort uvsorts vtrials vsorts labels sort validtrials tmpv tmpv2 subspace r5pc plotmin plotmax kmslast tmpunsortedidx tmpsortidx
set(handles.status, 'String', 'Wait'); drawnow
kms = str2num(get(handles.numcluster, 'String'));

if kms ~= kmslast
    tmpsortidx = zeros(1, size(tmpv2,2));
    tmpunsortedidx = find(tmpsortidx==0);
    kmslast = kms;
end

if get(handles.randomforrest, 'Value') == 1
    set(handles.kmeanscluster, 'Value', 1)
    set(handles.spectralcluster, 'Value',0)
    set(handles.kmedoidscluster, 'Value',0)
    set(handles.randomforrest, 'Value',0)
end

subspace = 2;
% tmpv = cat(2,vv{uvsorts, ch, :});
if kms == 1
    labels = ones(size(tmpv2(:,tmpunsortedidx),2),1);
%     labels = ones(size(tmpv2,2),1);
else

%     covtv = tmpv2*tmpv2';
%     [U,S,V]= svds(covtv,str2num(get(handles.numpcs, 'String')));
% 
%     ktmpv = (U'*tmpv2);
    ktmpv = (r5pc{ch}(:,1:str2num(get(handles.numpcs, 'String')))'*tmpv2(:,tmpunsortedidx));
%     ktmpv = (r5pc{ch}(:,1:str2num(get(handles.numpcs, 'String')))'*tmpv2);
    
    labels = cluster(ktmpv', kms, handles);
end

for ll = unique(labels)'
    tmpsortidx(tmpunsortedidx(labels==ll)) = ll;
end

for kk = 1:6
    hh = ['handles.axes' num2str(kk)];
    bdf = get(eval(hh), 'ButtonDownFcn');
    cla(eval(hh), 'reset')
    plot(eval(hh), tmpv2(:,tmpunsortedidx(labels==kk))); hold(eval(hh)); plot(eval(hh),mean(tmpv2(:,tmpunsortedidx(labels==kk)),2), 'k', 'linewidth', 4);
%     plot(eval(hh), tmpv2(:,(labels==kk))); hold(eval(hh)); plot(eval(hh),mean(tmpv2(:,(labels==kk)),2), 'k', 'linewidth', 4);
    set(eval(hh), 'YLim', [plotmin plotmax]);
    set(eval(hh), 'ButtonDownFcn', bdf)
end

if strcmp(PF.extradata{2}.data, 'fp30.index')
    scorefunc(handles)
end
set(handles.status, 'String', 'Ready'); drawnow
% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv tmpv2 subspace tmpclusters sortparams ens
global r5pc s5pc vs tridx tmps2 tmpsortidx tmpunsortedidx retry sc
set(handles.status, 'String', 'Wait'); drawnow

kms = str2num(get(handles.numcluster, 'String'));


% for ll = unique(labels)'
%     tmpsortidx(tmpunsortedidx(labels==ll)) = ll;
% end

for vt = validtrials
    for ll = unique(labels)'
%         if get(eval(['handles.singleunit' num2str(ll)]), 'Value')
%             PF.rec(vt).v_SU{ll,ch} = tmpv2(:,((tridx== vt) &  (tmpsortidx' == ll)));
%             PF.rec(vt).plexon_spikes_SU{ll,ch} = tmps2(((tridx== vt) &  (tmpsortidx' == ll)));
%             PF.rec(vt).F_SU{ll,ch} = sc(ll);
%         else
%             PF.rec(vt).v_MU{ll,ch} = tmpv2(:,((tridx== vt) &  (tmpsortidx' == ll)));
%             PF.rec(vt).plexon_spikes_MU{ll,ch} = tmps2(((tridx== vt) &  (tmpsortidx' == ll)));
%             PF.rec(vt).F_MU{ll,ch} = sc(ll);
%         end
        PF.rec(vt).v_final{ll,ch} = tmpv2(:,((tridx== vt) &  (tmpsortidx' == ll)));
        PF.rec(vt).plexon_spikes_final{ll,ch} = tmps2(((tridx== vt) &  (tmpsortidx' == ll)));
        PF.rec(vt).F_final{ll,ch} = sc(ll);
        if get(eval(['handles.singleunit' num2str(ll)]), 'Value')
            PF.rec(vt).type_final{ll,ch} = 'single';
        elseif get(eval(['handles.multiunit' num2str(ll)]), 'Value')
            PF.rec(vt).type_final{ll,ch} = 'multi';
        else
            PF.rec(vt).type_final{ll,ch} = 'noise';
        end
    end

    if size(PF.rec(vt).v_final(:,ch),1) > ll
        PF.rec(vt).v_final{ll+1,ch} = [];
        PF.rec(vt).plexon_spikes_final{ll+1,ch} = [];
        PF.rec(vt).F_final{ll+1,ch} = [];
        PF.rec(vt).type_final{ll+1,ch} = [];
    end

    if size(PF.rec(vt).v_final(:,ch),1) > ll+1
        PF.rec(vt).v_final{ll+2,ch} = [];
        PF.rec(vt).plexon_spikes_final{ll+2,ch} = [];
        PF.rec(vt).F_final{ll+2,ch} = [];
        PF.rec(vt).type_final{ll+2,ch} = [];
    end

    if size(PF.rec(vt).v_final(:,ch),1) > ll+2
        PF.rec(vt).v_final{ll+3,ch} = [];
        PF.rec(vt).plexon_spikes_final{ll+3,ch} = [];
        PF.rec(vt).F_final{ll+3,ch} = [];
        PF.rec(vt).type_final{ll+3,ch} = [];
    end
end


if isempty(ens) || isempty(ens{ch})
    if kms>1
        %% need to make sure this is not keeping original category
%         tcv = cat(3, tmpclusters(validtrials).v);
%         tcv = tcv(:,ch,:);
        cnt=0;
        n = sum(tmpsortidx>0);
        trnData = zeros(30+29+30+10+2,n);
        Y = zeros(1,n);
        for ii = 1:kms

%             trn = cat(2,tcv{ii,:});
            trn = tmpv2(:,tmpsortidx==ii);
            idx=cnt+1:cnt+size(trn,2);
            trnData(1:30,idx) = trn;
            trnData(31:59,idx) = diff(trn);
            trnData(60:89,idx) = abs(fft(trn));
            trnData(90:94,idx) = r5pc{ch}(:,1:5)'*trn;
            trnData(95:99,idx) = s5pc{ch}(:,1:5)'*trnData(31:59,idx);
            heights = max(trn) - min(trn);
            trnData(100,idx) = heights;
            [minr,minc] = find(bsxfun(@(x,y) (x == y), trn, min(trn))==1);
            [maxr,maxc] = find(bsxfun(@(x,y) (x == y), trn, max(trn))==1);
            minr = minr(unique(minc));
            maxr = maxr(unique(maxc));
            widths = minr - maxr;
            trnData(101,idx) = widths;
            Y(idx) = ii;
            cnt = cnt+size(trn,2);
        end

        ens{ch} = TreeBagger(300, trnData', Y');
    else
        ens{ch}=1;
    end

end
retry=0;
delete(handles.figure1);




% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes on button press in peaks.
function peaks_Callback(hObject, eventdata, handles)
% hObject    handle to peaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global PF vv ch currsort uvsorts vtrials vsorts labels sort validtrials tmpv tmpv2 subspace r5pc plotmin plotmax kmslast tmpsortidx tmpunsortedidx
set(handles.status, 'String', 'Wait'); drawnow
kms = str2num(get(handles.numcluster, 'String'));

if get(handles.randomforrest, 'Value') == 1
    set(handles.kmeanscluster, 'Value', 1)
    set(handles.spectralcluster, 'Value',0)
    set(handles.kmedoidscluster, 'Value',0)
    set(handles.randomforrest, 'Value',0)
end


if kms ~= kmslast
    tmpsortidx = zeros(1, size(tmpv2,2));
    tmpunsortedidx = find(tmpsortidx==0);
    kmslast = kms;
end


subspace = 1;


if kms == 1
    labels = ones(size(tmpv2,2),1);
else

    heights = max(tmpv2(:,tmpunsortedidx)) - min(tmpv2(:,tmpunsortedidx));
    heights = heights - mean(heights);
    heights = heights ./ std(heights);
    [minr,minc] = find(bsxfun(@(x,y) (x == y), tmpv2(:,tmpunsortedidx), min(tmpv2(:,tmpunsortedidx)))==1);
    [maxr,maxc] = find(bsxfun(@(x,y) (x == y), tmpv2(:,tmpunsortedidx), max(tmpv2(:,tmpunsortedidx)))==1);
    minr = minr(unique(minc));
    maxr = maxr(unique(maxc));
    widths = minr - maxr;
    widths = widths - mean(widths);
    widths = widths ./ std(widths);

    % 
%     covtv = tmpv2*tmpv2';
%     [U,S,V]= svds(covtv,str2num(get(handles.numpcs, 'String')));

    ktmpv2 = (r5pc{ch}(:,1:str2num(get(handles.numpcs, 'String')))'*tmpv2(:,tmpunsortedidx));

    ktmpv = [widths'; heights; ktmpv2];


    labels = cluster(ktmpv', kms, handles);
end

for ll = unique(labels)'
    tmpsortidx(tmpunsortedidx(labels==ll)) = ll;
end

for kk = 1:6
    hh = ['handles.axes' num2str(kk)];
    bdf = get(eval(hh), 'ButtonDownFcn');
    cla(eval(hh), 'reset')
    plot(eval(hh), tmpv2(:,(labels==kk))); hold(eval(hh)); plot(eval(hh),mean(tmpv2(:,(labels==kk)),2), 'k', 'linewidth', 4);
    set(eval(hh), 'YLim', [plotmin plotmax]);
    set(eval(hh), 'ButtonDownFcn', bdf)
end
if strcmp(PF.extradata{2}.data, 'fp30.index')
    scorefunc(handles)
end
set(handles.status, 'String', 'Ready'); drawnow


function labels = cluster(data, k, handles)
global ch tmpv2 validtrials subspace sortparams
global heightsmean heightsstd widthsmean widthsstd ens r5pc s5pc kms tmpsortidx tmpunsortedidx clusize vs tridx tmps2

if k ==1 | k ==0
    labels = ones(size(data,1),1);
    ens{ch} = 1;
else

    if get(handles.spectralcluster, 'Value') == 1
    % 
        ens{ch} = [];
        W = SimGraph(data', 1);
        A = SpectralClustering(W, k, 3);
        labels = full(sum(bsxfun(@times, A, 1:kms),2));
    elseif get(handles.kmeanscluster, 'Value') == 1
        ens{ch} = [];
        if size(data,1)<size(data,2)
            labels = fkmeans(data,k);
        else
            try

                centers = zeros(k,size(data,2));
                for ii = 1:k

                    trn = tmpv2(:,tmpsortidx==ii);

                    switch subspace
                        case 1

                            heights = max(trn) - min(trn);
                            heights = heights - heightsmean;
                            heights = heights ./ heightsstd;
                            [minr,minc] = find(bsxfun(@(x,y) (x == y), trn, min(trn))==1);
                            [maxr,maxc] = find(bsxfun(@(x,y) (x == y), trn, max(trn))==1);
                            minr = minr(unique(minc));
                            maxr = maxr(unique(maxc));
                            widths = minr - maxr;
                            widths = widths - widthsmean;
                            widths = widths ./ widthsstd;

                            trn2 = (r5pc{ch}(:,1:str2num(get(handles.numpcs, 'String')))'*trn);

                            trn = [widths'; heights; trn2];
                        case 2

                            trn = r5pc{ch}(:,1:str2num(get(handles.numpcs, 'String')))'*trn;
                    end

                    centers(ii,:)= mean(trn,2)';
                end

                labels = kmeans(data, k, 'start', centers);
            catch
                labels = kmeans(data, k);
            end
            
        end
    elseif get(handles.kmedoidscluster, 'Value') == 1
        ens{ch} = [];
        if size(data,1)==1
            labels=1;
        else
            labels = kmedoids(data', k)';
        end
    elseif get(handles.randomforrest, 'Value') == 1
        
        if exist('sortparams') && ~isempty(sortparams) && ~isempty(sortparams.ens{ch})
            X = zeros(101, size(data,2));
            X(1:30,:) = data;
            X(31:59,:) = diff(data);
            X(60:89,:) = abs(fft(data));
            X(90:94,:) = sortparams.r5pc{ch}(:,1:5)'*data;
            X(95:99,:) = sortparams.s5pc{ch}(:,1:5)'*X(31:59,:);
            heights = max(data) - min(data);
            X(100,:) = heights;
            [minr,minc] = find(bsxfun(@(x,y) (x == y), data, min(data))==1);
            [maxr,maxc] = find(bsxfun(@(x,y) (x == y), data, max(data))==1);
            minr = minr(unique(minc));
            maxr = maxr(unique(maxc));
            widths = minr - maxr;
            X(101,:) = widths;

            labels = predict(sortparams.ens{ch}, X');
            labels = str2num(cell2mat(labels));
        else
            try
                
%                 tcv = cat(3, tmpclusters(validtrials).v);
%                 tcv = tcv(:,ch,:);
                cnt=0;
                n = sum(tmpsortidx>0);
                trnData = zeros(30+29+30+10+2,n);
                Y = zeros(1,n);
                for ii = 1:k

%                     trn = cat(2,tcv{ii,:});
                    trn = tmpv2(:,tmpsortidx==ii);
                    idx=cnt+1:cnt+size(trn,2);
                    trnData(1:30,idx) = trn;
                    trnData(31:59,idx) = diff(trn);
                    trnData(60:89,idx) = abs(fft(trn));
                    trnData(90:94,idx) = r5pc{ch}(:,1:5)'*trn;
                    trnData(95:99,idx) = s5pc{ch}(:,1:5)'*trnData(31:59,idx);
                    heights = max(trn) - min(trn);
                    trnData(100,idx) = heights;
                    [minr,minc] = find(bsxfun(@(x,y) (x == y), trn, min(trn))==1);
                    [maxr,maxc] = find(bsxfun(@(x,y) (x == y), trn, max(trn))==1);
                    minr = minr(unique(minc));
                    maxr = maxr(unique(maxc));
                    widths = minr - maxr;
                    trnData(101,idx) = widths;
                    Y(idx) = ii;
                    cnt = cnt+size(trn,2);
                end

                ens{ch} = TreeBagger(300, trnData', Y');
                if ~isempty(data)
%                     data = data';

                    X = zeros(101, size(data,2));
                    X(1:30,:) = data;
                    X(31:59,:) = diff(data);
                    X(60:89,:) = abs(fft(data));
                    X(90:94,:) = r5pc{ch}(:,1:5)'*data;
                    X(95:99,:) = s5pc{ch}(:,1:5)'*X(31:59,:);
                    heights = max(data) - min(data);
                    X(100,:) = heights;
                    [minr,minc] = find(bsxfun(@(x,y) (x == y), data, min(data))==1);
                    [maxr,maxc] = find(bsxfun(@(x,y) (x == y), data, max(data))==1);
                    minr = minr(unique(minc));
                    maxr = maxr(unique(maxc));
                    widths = minr - maxr;
                    X(101,:) = widths;

                    labels = predict(ens{ch}, X');
                    labels = str2num(cell2mat(labels));
                else
                    labels = [];
                end
            catch
                
                error('There is no training data for the RF')
                
            end
        end

    end
end

% --- Executes on button press in kmeanscluster.
function kmeanscluster_Callback(hObject, eventdata, handles)
% hObject    handle to kmeanscluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of kmeanscluster
set(hObject, 'Value', 1)
set(handles.spectralcluster, 'Value',0)
set(handles.kmedoidscluster, 'Value',0)
set(handles.randomforrest, 'Value',0)

% --- Executes on button press in spectralcluster.
function spectralcluster_Callback(hObject, eventdata, handles)
% hObject    handle to spectralcluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of spectralcluster
set(hObject, 'Value', 1)
set(handles.kmeanscluster, 'Value',0)
set(handles.kmedoidscluster, 'Value',0)
set(handles.randomforrest, 'Value',0)

% --- Executes on button press in kmedoidscluster.
function kmedoidscluster_Callback(hObject, eventdata, handles)
% hObject    handle to kmedoidscluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of kmedoidscluster
set(hObject, 'Value', 1)
set(handles.spectralcluster, 'Value',0)
set(handles.kmeanscluster, 'Value',0)
set(handles.randomforrest, 'Value',0)


% --- Executes on button press in randomforrest.
function randomforrest_Callback(hObject, eventdata, handles)
% hObject    handle to randomforrest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of randomforrest
set(hObject, 'Value', 1)
set(handles.spectralcluster, 'Value',0)
set(handles.kmedoidscluster, 'Value',0)
set(handles.kmeanscluster, 'Value',0)


% --- Executes on button press in score.
function score_Callback(hObject, eventdata, handles)
% hObject    handle to score (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

scorefunc(handles)

function scorefunc(handles)
global PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv tmpv2 subspace tmpclusters sortparams tmpsortidx tmpunsortedidx clusize vs tridx tmps2 sc
kms = str2num(get(handles.numcluster, 'String'));
set(handles.status, 'String', 'Wait'); drawnow
tmpclusters2=tmpclusters;
PF2 = PF;
labels2 =labels;
%%%
tmpsortidx2 = tmpsortidx;

if isempty(labels)
    labels = tmpsortidx';
else

    for ll = unique(labels)'
        tmpsortidx2(tmpunsortedidx(labels==ll)) = ll;
    end
end

for vt = validtrials
    for ll = unique(labels)'
        PF2.rec(vt).v_final{ll} = tmpv2(:,((tridx== vt) &  (tmpsortidx2' == ll)));
        PF2.rec(vt).plexon_spikes_final{ll} = tmps2(((tridx== vt) &  (tmpsortidx2' == ll)));
    end
end



[rasters, FrameRate_Hz] = tmpraster(PF2);


combine = [];
sc = [];
for kk = 1:6
    
    hh = ['handles.score' num2str(kk)];
    hh2 = ['handles.checkbox' num2str(kk)];
    hh3 = ['handles.raster' num2str(kk)];
    if kk <= kms
        %sc = avgCorrMetric(rasters{kk,ch});
        if length(rasters)>=kk
            nansum(nansum(rasters{kk}))
            r1 = compress_raster(rasters{kk});
            %sc=calcExplainableVariance2(r1);
            sc(kk)=calcFstatistics2(r1);
            imagesc(r1','parent',eval(hh3)); colormap gray
            set(eval(hh), 'String', ['F: ' num2str(sc(kk))])

            if get(eval(hh2), 'Value') == 1
                combine = [combine kk];
            end
        end
    else
        imagesc([],'parent',eval(hh3)); colormap gray;
        set(eval(hh), 'String', [])
    end
end

rr=0;

for kk = combine
    rr = rr + compress_raster(rasters{kk});
end

if ~isempty(rr)
    %sc=calcExplainableVariance2(rr);
    scc=calcFstatistics2(rr);
    

    set(handles.combined, 'String', ['Combined F: ' num2str(scc)]);
else
    set(handles.combined, 'String', []);
end

set(handles.status, 'String', 'Ready'); drawnow


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6


% --- Executes on button press in discard.
function discard_Callback(hObject, eventdata, handles)
% hObject    handle to discard (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ens ch retry
ens{ch} = [];
retry=0;
delete(handles.figure1);


% --- Executes on button press in pcplot.
function pcplot_Callback(hObject, eventdata, handles)
% hObject    handle to pcplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ch tmpclusters validtrials subspace sortparams labels
global heightsmean heightsstd widthsmean widthsstd ens r5pc s5pc kms tmpv2
figure;
colors = {'r', 'g', 'b', 'y', 'c', 'k'}

for kk = 1:kms
    tmp = r5pc{ch}(:,1:2)'*tmpv2(:,labels==kk);
    hold on; scatter(tmp(1,:), tmp(2,:), colors{kk})
end


% --- Executes on button press in retry.
function retry_Callback(hObject, eventdata, handles)
% hObject    handle to retry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global retry ch ens
ens{ch} = [];
retry=1;
delete(handles.figure1);

function [plotmin, plotmax] = findplotminmax(tmpv2)

plotokflag = 0;
cnt = 0;

plotmin = min(tmpv2(:));
plotmax = max(tmpv2(:));

while plotokflag == 0
    

    [~,maxidx] = find(tmpv2 == max(tmpv2(:)));
    [~,minidx] = find(tmpv2 == min(tmpv2(:)));

    if plotmax >  max(4*std(tmpv2'))
        tmpv2(:,maxidx) = [];
        plotmin = min(tmpv2(:));
        plotmax = max(tmpv2(:));
    else
        plotokflag=1;
    end
    
    cnt = cnt+1;
    if cnt>10
        plotokfflag=1;
    end
end


% --- Executes on button press in multiunit1.
function multiunit1_Callback(hObject, eventdata, handles)
% hObject    handle to multiunit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of multiunit1
set(handles.singleunit1, 'Value',0)
set(handles.noise1, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in singleunit1.
function singleunit1_Callback(hObject, eventdata, handles)
% hObject    handle to singleunit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of singleunit1
set(handles.multiunit1, 'Value',0)
set(handles.noise1, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in multiunit2.
function multiunit2_Callback(hObject, eventdata, handles)
% hObject    handle to multiunit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of multiunit2
set(handles.singleunit2, 'Value',0)
set(handles.noise2, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in singleunit2.
function singleunit2_Callback(hObject, eventdata, handles)
% hObject    handle to singleunit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of singleunit2
set(handles.multiunit2, 'Value',0)
set(handles.noise2, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in multiunit3.
function multiunit3_Callback(hObject, eventdata, handles)
% hObject    handle to multiunit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of multiunit3
set(handles.singleunit3, 'Value',0)
set(handles.noise3, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in singleunit3.
function singleunit3_Callback(hObject, eventdata, handles)
% hObject    handle to singleunit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of singleunit3
set(handles.multiunit3, 'Value',0)
set(handles.noise3, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in multiunit4.
function multiunit4_Callback(hObject, eventdata, handles)
% hObject    handle to multiunit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of multiunit4
set(handles.singleunit4, 'Value',0)
set(handles.noise4, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in singleunit4.
function singleunit4_Callback(hObject, eventdata, handles)
% hObject    handle to singleunit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of singleunit4
set(handles.multiunit4, 'Value',0)
set(handles.noise4, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in multiunit5.
function multiunit5_Callback(hObject, eventdata, handles)
% hObject    handle to multiunit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of multiunit5
set(handles.singleunit5, 'Value',0)
set(handles.noise5, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in singleunit5.
function singleunit5_Callback(hObject, eventdata, handles)
% hObject    handle to singleunit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of singleunit5
set(handles.multiunit5, 'Value',0)
set(handles.noise5, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in multiunit6.
function multiunit6_Callback(hObject, eventdata, handles)
% hObject    handle to multiunit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of multiunit6
set(handles.singleunit6, 'Value',0)
set(handles.noise6, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in singleunit6.
function singleunit6_Callback(hObject, eventdata, handles)
% hObject    handle to singleunit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of singleunit6
set(handles.multiunit6, 'Value',0)
set(handles.noise6, 'Value',0)
set(hObject, 'Value',1)


% --- Executes on button press in noise1.
function noise1_Callback(hObject, eventdata, handles)
% hObject    handle to noise1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of noise1
set(handles.singleunit1, 'Value',0)
set(handles.multiunit1, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in noise2.
function noise2_Callback(hObject, eventdata, handles)
% hObject    handle to noise2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of noise2
set(handles.singleunit2, 'Value',0)
set(handles.multiunit2, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in noise3.
function noise3_Callback(hObject, eventdata, handles)
% hObject    handle to noise3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of noise3
set(handles.singleunit3, 'Value',0)
set(handles.multiunit3, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in noise4.
function noise4_Callback(hObject, eventdata, handles)
% hObject    handle to noise4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of noise4
set(handles.singleunit4, 'Value',0)
set(handles.multiunit4, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in noise5.
function noise5_Callback(hObject, eventdata, handles)
% hObject    handle to noise5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of noise5
set(handles.singleunit5, 'Value',0)
set(handles.multiunit5, 'Value',0)
set(hObject, 'Value',1)

% --- Executes on button press in noise6.
function noise6_Callback(hObject, eventdata, handles)
% hObject    handle to noise6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of noise6
set(handles.singleunit6, 'Value',0)
set(handles.multiunit6, 'Value',0)
set(hObject, 'Value',1)
