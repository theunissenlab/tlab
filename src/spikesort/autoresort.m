function varargout = autoresort(varargin)
% AUTORESORT MATLAB code for autoresort.fig
%      AUTORESORT, by itself, creates a new AUTORESORT or raises the existing
%      singleton*.
%
%      H = AUTORESORT returns the handle to a new AUTORESORT or the handle to
%      the existing singleton*.
%
%      AUTORESORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AUTORESORT.M with the given input arguments.
%
%      AUTORESORT('Property','Value',...) creates a new AUTORESORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before autoresort_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to autoresort_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help autoresort

% Last Modified by GUIDE v2.5 27-Mar-2012 17:58:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @autoresort_OpeningFcn, ...
                   'gui_OutputFcn',  @autoresort_OutputFcn, ...
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


% --- Executes just before autoresort is made visible.
function autoresort_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to autoresort (see VARARGIN)
global PF vv ch currsort uvsorts vtrials vsorts labels sort validtrials tmpv tmpv2 tmpclusters sortparams kms plotmin plotmax
% Choose default command line output for autoresort
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes autoresort wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% kms = str2num(get(handles.numcluster, 'String'));


        
% --- Outputs from this function are returned to the command line.
function varargout = autoresort_OutputFcn(hObject, eventdata, handles) 
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


% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv tmpv2 subspace tmpclusters sortparams ens
global r5pc s5pc sortidx
set(handles.status, 'String', 'Wait'); drawnow
% trsz = cat(1,squeeze(cellfun(@(x) size(x,2), vv(uvsorts, ch, :))));
% trsz = sum(trsz);
% 
% % tmpv2 = cat(2,vv{uvsorts, ch, :});
% pv = cat(3,PF.rec(:).plexon_spikes);
% tmppv2 = cat(2,pv{uvsorts, ch, :});
% cnt=0;
% for vt = unique(vtrials)'
%     cnt = cnt+1;
%     for ll = unique(labels)'
%         [tvidx] = find(labels==ll);
%         svidx = tvidx(tvidx <= trsz(cnt));
%         
%         PF.rec(validtrials(vt)).v_final{ll,ch} = tmpv2(:,svidx);
%         PF.rec(validtrials(vt)).plexon_spikes_final{ll,ch} = tmppv2(svidx);
%     end
%     
%     labels(1:trsz(cnt)) = [];
%     tmpv2(:,1:trsz(cnt)) = [];
%     tmppv2(1:trsz(cnt)) = [];
%     
% end

kms = str2num(get(handles.numcluster, 'String'));

%%%
for ll = unique(labels)'
    [tvidx] = find(labels==ll);
    labels(labels==ll) = [];
    
    for ccsort = uvsorts'
        
        vv = cat(3,PF.rec(:).v);

        trsz = cat(1,squeeze(cellfun(@(x) size(x,2), vv(ccsort, ch, vtrials))));
        trsz(vsorts~=ccsort)=0;

        cnt = 0;

        for vt = vtrials'
            cnt = cnt+1;
            delidx = tvidx(tvidx <= trsz(cnt));
            tvidx = tvidx(tvidx > trsz(cnt));
            tvidx = tvidx - trsz(cnt);
            if size(tmpclusters(validtrials(vt)).v,1) < ll
                tmpclusters(validtrials(vt)).v{ll,ch} = [];
                tmpclusters(validtrials(vt)).plexon_spikes{ll,ch} = [];
            end
            if size(tmpclusters(validtrials(vt)).v,2) < ch
                tmpclusters(validtrials(vt)).v{ll,ch} = [];
                tmpclusters(validtrials(vt)).plexon_spikes{ll,ch} = [];
            end

            tmpclusters(validtrials(vt)).v{ll,ch} = [tmpclusters(validtrials(vt)).v{ll,ch} PF.rec(validtrials(vt)).v{ccsort,ch}(:,delidx)];
            tmpclusters(validtrials(vt)).plexon_spikes{ll,ch} = [tmpclusters(validtrials(vt)).plexon_spikes{ll,ch} PF.rec(validtrials(vt)).plexon_spikes{ccsort,ch}(delidx)];
            [~,dsidx] = sort(tmpclusters(validtrials(vt)).plexon_spikes{ll,ch}, 'ascend');
            tmpclusters(validtrials(vt)).v{ll,ch} = tmpclusters(validtrials(vt)).v{ll,ch}(:,dsidx);
            tmpclusters(validtrials(vt)).plexon_spikes{ll,ch} = tmpclusters(validtrials(vt)).plexon_spikes{ll,ch}(dsidx);
            PF.rec(validtrials(vt)).v{ccsort,ch}(:,delidx) = [];
            PF.rec(validtrials(vt)).plexon_spikes{ccsort,ch}(delidx) = [];
        end
    end
end

if isempty(ens) || isempty(ens{ch})
    if kms>1
        
        tcv = cat(3, tmpclusters(validtrials).v);
        tcv = tcv(:,ch,:);
        cnt=0;
        n = size(cat(2,tcv{:}),2);
        trnData = zeros(30+29+30+10+2,n);
        Y = zeros(1,n);
        for ii = 1:kms

            trn = cat(2,tcv{ii,:});
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

sortidx = sortidx + 1;
makeplots(handles)



% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);




function labels = cluster(data, k, handles)
global ch tmpclusters validtrials subspace sortparams
global heightsmean heightsstd widthsmean widthsstd ens r5pc s5pc kms

if k ==1 | k ==0
    labels = ones(size(data,1),1);
    ens{ch} = 1;
else

    if get(handles.spectralcluster, 'Value') == 1
    % keyboard
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
                if exist('sortparams') && ~isempty(sortparams)
                    tcv = cat(3, sortparams.tmpclusters(validtrials).v);
                    tcv = tcv(:,ch,:);
                else
                    tcv = cat(3, tmpclusters(validtrials).v);
                    tcv = tcv(:,ch,:);
                end

                centers = zeros(k,size(data,2));
                for ii = 1:k
                    trn = cat(2,tcv{ii,:});

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
        
        if exist('sortparams') && ~isempty(sortparams)
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
                tcv = cat(3, tmpclusters(validtrials).v);
                tcv = tcv(:,ch,:);
                cnt=0;
                n = size(cat(2,tcv{:}),2);
                trnData = zeros(30+29+30+10+2,n);
                Y = zeros(1,n);
                for ii = 1:k

                    trn = cat(2,tcv{ii,:});
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

                data = data';

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
global PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv tmpv2 subspace tmpclusters sortparams
kms = str2num(get(handles.numcluster, 'String'));
set(handles.status, 'String', 'Wait'); drawnow
tmpclusters2=tmpclusters;
PF2 = PF;
labels2 =labels;
%%%
for ll = unique(labels)'
    [tvidx] = find(labels2==ll);
    labels2(labels2==ll) = [];
    
    for ccsort = uvsorts'
        
        vv = cat(3,PF2.rec(:).v);

        trsz = cat(1,squeeze(cellfun(@(x) size(x,2), vv(ccsort, ch, vtrials))));
        trsz(vsorts~=ccsort)=0;

        cnt = 0;

        for vt = vtrials'
            cnt = cnt+1;
            delidx = tvidx(tvidx <= trsz(cnt));
            tvidx = tvidx(tvidx > trsz(cnt));
            tvidx = tvidx - trsz(cnt);
            if size(tmpclusters2(validtrials(vt)).v,1) < ll
                tmpclusters2(validtrials(vt)).v{ll,ch} = [];
                tmpclusters2(validtrials(vt)).plexon_spikes{ll,ch} = [];
            end
            if size(tmpclusters2(validtrials(vt)).v,2) < ch
                tmpclusters2(validtrials(vt)).v{ll,ch} = [];
                tmpclusters2(validtrials(vt)).plexon_spikes{ll,ch} = [];
            end

            tmpclusters2(validtrials(vt)).v{ll,ch} = [tmpclusters2(validtrials(vt)).v{ll,ch} PF2.rec(validtrials(vt)).v{ccsort,ch}(:,delidx)];
            tmpclusters2(validtrials(vt)).plexon_spikes{ll,ch} = [tmpclusters2(validtrials(vt)).plexon_spikes{ll,ch} PF2.rec(validtrials(vt)).plexon_spikes{ccsort,ch}(delidx)];
            [~,dsidx] = sort(tmpclusters2(validtrials(vt)).plexon_spikes{ll,ch}, 'ascend');
            tmpclusters2(validtrials(vt)).v{ll,ch} = tmpclusters2(validtrials(vt)).v{ll,ch}(:,dsidx);
            tmpclusters2(validtrials(vt)).plexon_spikes{ll,ch} = tmpclusters2(validtrials(vt)).plexon_spikes{ll,ch}(dsidx);
            PF2.rec(validtrials(vt)).v{ccsort,ch}(:,delidx) = [];
            PF2.rec(validtrials(vt)).plexon_spikes{ccsort,ch}(delidx) = [];
        end
    end
end

for ii = 1:size(tmpclusters2,2)
    PF2.rec(ii).plexon_spikes = tmpclusters2(ii).plexon_spikes;
end

[rasters, FrameRate_Hz] = tmpraster(PF2);
combine = [];
for kk = 1:6
    
    hh = ['handles.score' num2str(kk)];
    hh2 = ['handles.checkbox' num2str(kk)];
    hh3 = ['handles.raster' num2str(kk)];
    if kk <= kms
        %sc = avgCorrMetric(rasters{kk,ch});
        nansum(nansum(rasters{kk,ch}))
        r1 = compress_raster(rasters{kk,ch});
        %sc=calcExplainableVariance2(r1);
        sc=calcFstatistics2(r1);
        imagesc(r1','parent',eval(hh3)); colormap gray
        set(eval(hh), 'String', ['EV: ' num2str(sc)])

        if get(eval(hh2), 'Value') == 1
            combine = [combine kk];
        end
    else
        imagesc([],'parent',eval(hh3)); colormap gray;
        set(eval(hh), 'String', [])
    end
end

rr=0;

for kk = combine
    rr = rr + compress_raster(rasters{kk,ch});
end

if ~isempty(rr)
    %sc=calcExplainableVariance2(rr);
    sc=calcFstatistics2(rr);
    

    set(handles.combined, 'String', ['Combined EV: ' num2str(sc)]);
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
global ens ch
ens{ch} = [];
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



function filename_Callback(hObject, eventdata, handles)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filename as text
%        str2double(get(hObject,'String')) returns contents of filename as a double


% --- Executes during object creation, after setting all properties.
function filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in autoresort.
function autoresort_Callback(hObject, eventdata, handles)
% hObject    handle to autoresort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.status, 'String', 'Wait'); drawnow
% set(handles.nextbutton, 'Enable', 'off')
global partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv goodflag tmpclusters fname sortparams
global heightsmean heightsstd widthsmean widthsstd s5pc r5pc ens sortedchans sortidx

fname = get(handles.filename, 'String');

if exist(fname) ~= 2
    error('File does not exist')
    return
end
    
load(fname)


for ii = 1:size(PF.rec,2)
    vsizes(ii) = size(PF.rec(ii).v,1);
end

fixidx = find(vsizes < max(vsizes));

for ii = fixidx
    if size(PF.rec(ii).v,1)>0
    
        PF.rec(ii).v{max(vsizes),1} = [];
    end
end

vv = cat(3,PF.rec(:).v);
vt = 1:size(PF.rec,2);
validtrials=[];
for ii = vt
    tmpclusters(ii).v{1,1}=[];
    tmpclusters(ii).plexon_spikes{1,1}=[];
    if isempty(PF.rec(ii).v)
        good = 0;
    else
        good = 1;
    end
    validtrials = [validtrials good];
end

validtrials = vt(logical(validtrials));


s1 = strfind(fname, '/');


sortparams = load([fname(1:s1(end)) 'sortparams.mat']);
ens = sortparams.ens;
s5pc = sortparams.s5pc; 
r5pc = sortparams.r5pc;
sortedchans = find(cellfun(@(x) ~isempty(x), sortparams.ens)==1);

sortidx=1;

makeplots(handles)

function makeplots(handles)
set(handles.status, 'String', 'Wait'); drawnow
% set(handles.nextbutton, 'Enable', 'off')
global partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv goodflag tmpclusters fname sortparams
global heightsmean heightsstd widthsmean widthsstd s5pc r5pc ens sortedchans sortidx
 
if sortidx > length(sortedchans)
    disp('Done')
    return
end

ch = sortedchans(sortidx);
set(handles.channelnum, 'String', ['Channel: ', num2str(ch)]); drawnow

[vsorts, vtrials] = find(cellfun(@(x) ~isempty(x),vv(:,ch,:))==1);

uvsorts = unique(vsorts);

tmpv2 = cat(2,vv{uvsorts, ch, :});
plotmin = min(tmpv2(:));
plotmax = max(tmpv2(:));

plotokflag = 0;

%         while plotokflag == 0
% 
%             [~,maxidx] = find(tmpv2 == max(tmpv2(:)));
%             [~,minidx] = find(tmpv2 == min(tmpv2(:)));
% 
%             if plotmax >  max(2*std(tmpv2'))
%                 tmpv2(:,maxidx) = [];
        plotmin = min(tmpv2(:));
        plotmax = max(tmpv2(:));
%             else
        plotokflag=1;
%             end
%         end

tmpv2 = cat(2,vv{uvsorts, ch, :});

if strcmp(class(sortparams.ens{ch}), 'double') && sortparams.ens{ch} == 1

    labels = ones(size(tmpv2,2),1);
    kms=1;
else
    X = zeros(101, size(tmpv2,2));
    X(1:30,:) = tmpv2;
    X(31:59,:) = diff(tmpv2);
    X(60:89,:) = abs(fft(tmpv2));
    X(90:94,:) = sortparams.r5pc{ch}(:,1:5)'*tmpv2;
    X(95:99,:) = sortparams.s5pc{ch}(:,1:5)'*X(31:59,:);
    heights = max(tmpv2) - min(tmpv2);
    X(100,:) = heights;
    [minr,minc] = find(bsxfun(@(x,y) (x == y), tmpv2, min(tmpv2))==1);
    [maxr,maxc] = find(bsxfun(@(x,y) (x == y), tmpv2, max(tmpv2))==1);
    minr = minr(unique(minc));
    maxr = maxr(unique(maxc));
    widths = minr - maxr;
    X(101,:) = widths;

    labels = predict(sortparams.ens{ch}, X');
    labels = str2num(cell2mat(labels));

    kms = str2num(sortparams.ens{ch}.ClassNames{end});
end

set(handles.numcluster, 'String', num2str(kms));

for kk = 1:6
    hh = ['handles.axes' num2str(kk)];
    bdf = get(eval(hh), 'ButtonDownFcn');
    cla(eval(hh), 'reset')
    plot(eval(hh), tmpv2(:,labels==kk)); set(eval(hh), 'YLim', [plotmin plotmax]);
    hold(eval(hh)); plot(eval(hh), mean(tmpv2(:,labels==kk),2), 'k', 'linewidth', 4)
    set(eval(hh), 'ButtonDownFcn', bdf)
end

set(handles.status, 'String', 'Ready'); drawnow


% --- Executes on button press in resort.
function resort_Callback(hObject, eventdata, handles)
% hObject    handle to resort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
