function varargout = resorter(varargin)
% RESORTER MATLAB code for resorter.fig
%      RESORTER, by itself, creates a new RESORTER or raises the existing
%      singleton*.
%
%      H = RESORTER returns the handle to a new RESORTER or the handle to
%      the existing singleton*.
%
%      RESORTER('CALLBACK',hObject,eventData,handles,...) calls the localff
%      function named CALLBACK in RESORTER.M with the given input arguments.
%
%      RESORTER('Property','Value',...) creates a new RESORTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before resorter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to resorter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help resorter

% Last Modified by GUIDE v2.5 25-Apr-2012 18:08:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @resorter_OpeningFcn, ...
                   'gui_OutputFcn',  @resorter_OutputFcn, ...
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


% --- Executes just before resorter is made visible.
function resorter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to resorter (see VARARGIN)

% Choose default command line output for resorter
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes resorter wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = resorter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



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


% --- Executes on button press in sortbutton.
function sortbutton_Callback(hObject, eventdata, handles)
% hObject    handle to sortbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
bdf2 = disable_all(handles);


clear PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv tmpv2 goodflag fname
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
        PF.rec(ii).plexon_spikes{max(vsizes),1}=[];
        PF.rec(ii).v{max(vsizes),1} = [];
    end
end

global partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv goodflag tmpclusters fname sortparams
global heightsmean heightsstd widthsmean widthsstd s5pc r5pc ens tmpsortidx clusize tmpv2 tmpunsortedidx vs tridx tmps2 retry undohistory colortable resorter_handles

resorter_handles = handles;

colortable{1} = [255 0 0];
colortable{2} = [255 128 0];
colortable{3} = [255 255 0];
colortable{4} = [0 255 0];
colortable{5} = [0 255 255];
colortable{6} = [0 128 255];
colortable{7} = [0 0 255];
colortable{8} = [128 0 255];
colortable{9} = [255 0 255];
colortable{10} = [255 0 128];
colortable{11} = [100 0 0];
colortable{12} = [100 100 0];
colortable{13} = [0 100 0];
colortable{14} = [0 100 100];
colortable{15} = [0 0 100];
colortable{16} = [100 0 100];

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
retry=0;


s1 = strfind(fname, '/');

if exist([fname(1:s1(end)) 'sortparams.mat'])
    sortparams = load([fname(1:s1(end)) 'sortparams.mat']);
    ens = sortparams.ens;
    r5pc = sortparams.r5pc;
    s5pc = sortparams.s5pc;
else
    sortparams = [];
end

    ch = 0;

if ~isempty(sortparams)
    finalcluster
else

    if isempty(currsort)
        currsort=0;
    end

    flag = 0;

    while flag == 0

        currsort = currsort + 1;

        if currsort > length(uvsorts)

            ch = ch + 1;
            currsort = 1;
            [vsorts, vtrials] = find(cellfun(@(x) ~isempty(x),vv(:,ch,:))==1);

            uvsorts = unique(vsorts);

        end

        kms = 6;

        set(handles.ChanNum, 'String', ['Channel: ' num2str(ch)]);


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
        
        
        
        
%         plotmin = min(tmpv2(:));
%         plotmax = max(tmpv2(:));
        tmpsortidx = zeros(1, size(tmpv2,2));
        tmpunsortedidx = find(tmpsortidx==0);
        
        [plotmin, plotmax] = findplotminmax(tmpv2);
        
        allclusters
        
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

%         if mean(tmpv(1,:))<mean(tmpv(7,:))
%             continue
%         else
            goodflag=1;
%         end

        if clusize(currsort) > str2num(get(handles.snipnum, 'string'));
            snpidx1 = 1;
            snpidx2 = str2num(get(handles.snipnum, 'string'));
            partialflag = 1;
        else
            snpidx1 = 1;
            snpidx2 = clusize(currsort);
            partialflag = 0;
        end
        
        undohistory{1}.tmpsortidx=tmpsortidx;
        undohistory{1}.snpidx1=snpidx1;
        undohistory{1}.snpidx2=snpidx2;
        undohistory{1}.currsort=currsort;
        undohistory{1}.partialflag = partialflag;
        snpidx = snpidx1:snpidx2;

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

% --- Executes on button press in nextbutton.
function nextbutton_Callback(hObject, eventdata, handles)
% hObject    handle to nextbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv tmpv2 goodflag r5pc s5pc
global heightsmean heightsstd widthsmean widthsstd ens tmpclusters fname tmpsortidx clusize tmpunsortedidx vs tridx tmps2 retry undohistory
bdf2 = disable_all(handles);
if isempty(ch)
    ch = 1;
end

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
%                 while isempty(intersect(snpidx1:snpidx2, tmpunsortedidx))
%                     currsort = currsort+1;
%                     snpidx1 = clusize(currsort-1)+1;
%                     if (clusize(currsort) - clusize(currsort-1)) > str2num(get(handles.snipnum, 'string'));
%                         snpidx2 = snpidx1 + str2num(get(handles.snipnum, 'string'))-1;
%                         partialflag = 1;
%                     else
%                         if currsort == length(uvsorts)
%                             snpidx2 = size(tmpv2,2);
%                         else
%                             snpidx2 = clusize(currsort);
%                         end
%                         partialflag = 0;
%                     end
%                 end
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
%         load(fname)
         s1 = strfind(fname, '/');
%         
%         for vt = validtrials
%             PF.rec(vt).v_final = tmpclusters(vt).v;
%             PF.rec(vt).plexon_spikes_final = tmpclusters(vt).plexon_spikes;
%         end
        
        save([fname(1:end-3) 'rs.mat'], 'PF')
        
        save([fname(1:s1(end)) 'sortparams.mat'], 'ens', 'heightsmean', 'heightsstd', 'widthsmean', 'widthsstd', 'r5pc', 's5pc')
        return
    end
            
    
    snpidx = intersect(snpidx1:snpidx2, tmpunsortedidx);

    
%     if mean(tmpv2(1,:))<mean(tmpv2(7,:))
% 
%         continue
%     else
        goodflag=1;
%     end

    if (currsort == 1 && snpidx1 == 1) || isempty(heightsmean)
%         keyboard
        [plotmin, plotmax] = findplotminmax(tmpv2);
        allclusters
        
%         tmpv2 = cat(2,vv{uvsorts, ch, :});
        

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

enable_all(handles, bdf2)
% --- Executes on button press in deletesc.
function deletesc_Callback(hObject, eventdata, handles)
% hObject    handle to deletesc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of deletesc
global clu
clu = NaN;
set(hObject, 'Value', 1)
set(handles.cluster2, 'Value',0)
set(handles.cluster3, 'Value',0)
set(handles.cluster4, 'Value',0)
set(handles.cluster5, 'Value',0)
set(handles.cluster6, 'Value',0)
set(handles.cluster1, 'Value',0)
set(handles.examine, 'Value',0)

% --- Executes on button press in examine.
function examine_Callback(hObject, eventdata, handles)
% hObject    handle to examine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of examine
set(hObject, 'Value', 1)
set(handles.cluster2, 'Value',0)
set(handles.cluster3, 'Value',0)
set(handles.cluster4, 'Value',0)
set(handles.cluster5, 'Value',0)
set(handles.cluster6, 'Value',0)
set(handles.deletesc, 'Value',0)
set(handles.cluster1, 'Value',0)

% --- Executes on mouse press over axes background.
function kmean1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to kmean1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotclick(1,handles)

% --- Executes on mouse press over axes background.
function kmean2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to kmean2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotclick(2,handles)
% --- Executes on mouse press over axes background.
function kmean3_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to kmean3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotclick(3,handles)

% --- Executes on mouse press over axes background.
function kmean4_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to kmean4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotclick(4,handles)

% --- Executes on mouse press over axes background.
function kmean5_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to kmean5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotclick(5,handles)

% --- Executes on mouse press over axes background.
function kmean6_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to kmean6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotclick(6,handles)


function pcnum_Callback(hObject, eventdata, handles)
% hObject    handle to pcnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pcnum as text
%        str2double(get(hObject,'String')) returns contents of pcnum as a double


% --- Executes during object creation, after setting all properties.
function pcnum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pcnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in peaks.
function kmeans_Callback(hObject, eventdata, handles)
% hObject    handle to peaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv2 vv
bdf2 = disable_all(handles);
% vv = cat(3,PF.rec(:).v);

kms = 6;

set(handles.ChanNum, 'String', ['Channel: ' num2str(ch)]);

% csort = uvsorts(currsort);

% tmpv = cat(2,vv{csort, ch, vtrials(vsorts==csort)});
% tmpv2 = cat(2,vv{uvsorts, ch, :});
% keyboard
plot(handles.mainplot, tmpv2(:,snpidx)); set(handles.mainplot, 'YLim', [plotmin plotmax]);
% keyboard

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
enable_all(handles, bdf2)


% --- Executes on button press in peaks.
function pca_Callback(hObject, eventdata, handles)
% hObject    handle to peaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global r5pc s5pc partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv2 vv

bdf2 = disable_all(handles);
% vv = cat(3,PF.rec(:).v);

kms = 6;

set(handles.ChanNum, 'String', ['Channel: ' num2str(ch)]);

% csort = uvsorts(currsort);
% 
% tmpv = cat(2,vv{csort, ch, vtrials(vsorts==csort)});
% tmpv2 = cat(2,vv{uvsorts, ch, :});

plot(handles.mainplot, tmpv2(:,snpidx)); set(handles.mainplot, 'YLim', [plotmin plotmax]);

if get(handles.slope, 'Value')==1

    if str2num(get(handles.pcnum, 'String'))>0
        ktmpv = (s5pc{ch}(:,str2num(get(handles.pcnum, 'String')))'*diff(tmpv2(:,snpidx)));
    else
        ktmpv = (s5pc{ch}(:,1)'*diff(tmpv2(:,snpidx)));
    end
else
    if str2num(get(handles.pcnum, 'String'))>0
        ktmpv = (r5pc{ch}(:,str2num(get(handles.pcnum, 'String')))'*tmpv2(:,snpidx));
    else
        ktmpv = (s5pc{ch}(:,1)'*diff(tmpv2(:,snpidx)));
    end
end

labels = cluster(ktmpv', kms, handles);


for kk = 1:kms
    hh = ['handles.kmean' num2str(kk)];
    bdf = get(eval(hh), 'ButtonDownFcn');
    cla(eval(hh), 'reset')
    plot(eval(hh), tmpv2(:,snpidx(labels==kk))); set(eval(hh), 'YLim', [plotmin plotmax]);
    set(eval(hh), 'ButtonDownFcn', bdf)
end
enable_all(handles, bdf2)


% --- Executes on button press in peaks.
function peaks_Callback(hObject, eventdata, handles)
% hObject    handle to peaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global r5pc s5pc partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv2
global heightsmean heightsstd widthsmean widthsstd
bdf2 = disable_all(handles);
% vv = cat(3,PF.rec(:).v);

kms = 6;

set(handles.ChanNum, 'String', ['Channel: ' num2str(ch)]);

% csort = uvsorts(currsort);

% tmpv = cat(2,vv{csort, ch, vtrials(vsorts==csort)});
% tmpv2 = cat(2,vv{uvsorts, ch, :});

plot(handles.mainplot, tmpv2(:,snpidx)); set(handles.mainplot, 'YLim', [plotmin plotmax]);

heights = max(tmpv2(:,snpidx)) - min(tmpv2(:,snpidx));
heights = heights - heightsmean;
heights = heights ./ heightsstd;
[minr,minc] = find(bsxfun(@(x,y) (x == y), tmpv2(:,snpidx), min(tmpv2(:,snpidx)))==1);
[maxr,maxc] = find(bsxfun(@(x,y) (x == y), tmpv2(:,snpidx), max(tmpv2(:,snpidx)))==1);
minr = minr(unique(minc));
maxr = maxr(unique(maxc));

widths = maxr - minr;
widths = widths - widthsmean;
widths = widths ./ widthsstd;

% keyboard

if get(handles.slope, 'Value')==1
    if str2num(get(handles.pcnum, 'String'))>0
        ktmpv2 = (s5pc{ch}(:,str2num(get(handles.pcnum, 'String')))'*diff(tmpv2(:,snpidx)));
    else
        ktmpv2 = [];
    end
else
    if str2num(get(handles.pcnum, 'String'))>0
        ktmpv2 = (r5pc{ch}(:,str2num(get(handles.pcnum, 'String')))'*tmpv2(:,snpidx));
    else
        ktmpv2 = [];
    end
end

ktmpv = [widths'; heights; ktmpv2];


labels = cluster(ktmpv', kms, handles);


for kk = 1:kms
    hh = ['handles.kmean' num2str(kk)];
    bdf = get(eval(hh), 'ButtonDownFcn');
    cla(eval(hh), 'reset')
    plot(eval(hh), tmpv2(:,snpidx(labels==kk))); set(eval(hh), 'YLim', [plotmin plotmax]);
    set(eval(hh), 'ButtonDownFcn', bdf)
end
enable_all(handles, bdf2)


% --- Executes on button press in kmeanscluster.
function kmeanscluster_Callback(hObject, eventdata, handles)
% hObject    handle to kmeanscluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of kmeanscluster
set(hObject, 'Value', 1)
set(handles.spectralcluster, 'Value',0)
set(handles.kmedoidscluster, 'Value',0)
    

% --- Executes on button press in spectralcluster.
function spectralcluster_Callback(hObject, eventdata, handles)
% hObject    handle to spectralcluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of spectralcluster
set(hObject, 'Value', 1)
set(handles.kmeanscluster, 'Value',0)
set(handles.kmedoidscluster, 'Value',0)

% --- Executes on button press in kmedoidscluster.
function kmedoidscluster_Callback(hObject, eventdata, handles)
% hObject    handle to kmedoidscluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of kmedoidscluster
set(hObject, 'Value', 1)
set(handles.spectralcluster, 'Value',0)
set(handles.kmeanscluster, 'Value',0)





% --- Executes on button press in raw.
function raw_Callback(hObject, eventdata, handles)
% hObject    handle to raw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of raw
set(hObject, 'Value', 1)
set(handles.slope, 'Value',0)


% --- Executes on button press in slope.
function slope_Callback(hObject, eventdata, handles)
% hObject    handle to slope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of slope
set(hObject, 'Value', 1)
set(handles.raw, 'Value',0)


% --- Executes on button press in cluster1.
function cluster1_Callback(hObject, eventdata, handles)
% hObject    handle to cluster1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cluster1
global clu
clu = 1;
set(hObject, 'Value', 1)
set(handles.cluster2, 'Value',0)
set(handles.cluster3, 'Value',0)
set(handles.cluster4, 'Value',0)
set(handles.cluster5, 'Value',0)
set(handles.cluster6, 'Value',0)
set(handles.deletesc, 'Value',0)
set(handles.examine, 'Value',0)


% --- Executes on button press in cluster2.
function cluster2_Callback(hObject, eventdata, handles)
% hObject    handle to cluster2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cluster2
global clu
clu = 2;
set(hObject, 'Value', 1)
set(handles.cluster1, 'Value',0)
set(handles.cluster3, 'Value',0)
set(handles.cluster4, 'Value',0)
set(handles.cluster5, 'Value',0)
set(handles.cluster6, 'Value',0)
set(handles.deletesc, 'Value',0)
set(handles.examine, 'Value',0)

% --- Executes on button press in cluster3.
function cluster3_Callback(hObject, eventdata, handles)
% hObject    handle to cluster3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cluster3
global clu
clu = 3;
set(hObject, 'Value', 1)
set(handles.cluster2, 'Value',0)
set(handles.cluster1, 'Value',0)
set(handles.cluster4, 'Value',0)
set(handles.cluster5, 'Value',0)
set(handles.cluster6, 'Value',0)
set(handles.deletesc, 'Value',0)
set(handles.examine, 'Value',0)

% --- Executes on button press in cluster4.
function cluster4_Callback(hObject, eventdata, handles)
% hObject    handle to cluster4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cluster4
global clu
clu = 4;
set(hObject, 'Value', 1)
set(handles.cluster2, 'Value',0)
set(handles.cluster3, 'Value',0)
set(handles.cluster1, 'Value',0)
set(handles.cluster5, 'Value',0)
set(handles.cluster6, 'Value',0)
set(handles.deletesc, 'Value',0)
set(handles.examine, 'Value',0)

% --- Executes on button press in cluster5.
function cluster5_Callback(hObject, eventdata, handles)
% hObject    handle to cluster5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cluster5
global clu
clu = 5;
set(hObject, 'Value', 1)
set(handles.cluster2, 'Value',0)
set(handles.cluster3, 'Value',0)
set(handles.cluster4, 'Value',0)
set(handles.cluster1, 'Value',0)
set(handles.cluster6, 'Value',0)
set(handles.deletesc, 'Value',0)
set(handles.examine, 'Value',0)

% --- Executes on button press in cluster6.
function cluster6_Callback(hObject, eventdata, handles)
% hObject    handle to cluster6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cluster6
global clu
clu = 6;
set(hObject, 'Value', 1)
set(handles.cluster2, 'Value',0)
set(handles.cluster3, 'Value',0)
set(handles.cluster4, 'Value',0)
set(handles.cluster5, 'Value',0)
set(handles.cluster1, 'Value',0)
set(handles.deletesc, 'Value',0)
set(handles.examine, 'Value',0)


function snipnum_Callback(hObject, eventdata, handles)
% hObject    handle to snipnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of snipnum as text
%        str2double(get(hObject,'String')) returns contents of snipnum as a double


% --- Executes during object creation, after setting all properties.
function snipnum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to snipnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plotclick(labelnum, handles)
global partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels validtrials tmpv tmpv2 tmpclusters r5pc s5pc 
global heightsmean heightsstd widthsmean widthsstd goodflag fname tmpsortidx tmpunsortedidx clusize vs tridx tmps2 retry undohistory clu

bdf2 = disable_all(handles);
% if get(handles.deletesc, 'Value') == 1
% 
%     tmpsortidx(snpidx(find(labels==labelnum))) = NaN;
%     tmpunsortedidx = find(tmpsortidx==0);
% 
% 
%     clu = 0;
% elseif get(handles.cluster1, 'Value') == 1
%     
%     clu=1;
% elseif get(handles.cluster2, 'Value') == 1
%     
%     clu=2;
% elseif get(handles.cluster3, 'Value') == 1
%     
%     clu=3;
% elseif get(handles.cluster4, 'Value') == 1
%     
%     clu=4;
% elseif get(handles.cluster5, 'Value') == 1
%     
%     clu=5;
% elseif get(handles.cluster6, 'Value') == 1
%     
%     clu=6;
% else
%     clu=0;
% end

tmpsortidx(snpidx(find(labels==labelnum))) = clu;
tmpunsortedidx = find(tmpsortidx==0);

if clu>0

    tempcluster
    set(handles.figure1, 'Visible', 'on')
end

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
                
%                 while isempty(intersect(snpidx1:snpidx2, tmpunsortedidx))
%                     currsort = currsort+1;
%                     snpidx1 = clusize(currsort-1)+1;
%                     if (clusize(currsort) - clusize(currsort-1)) > str2num(get(handles.snipnum, 'string'));
%                         snpidx2 = snpidx1 + str2num(get(handles.snipnum, 'string'))-1;
%                         partialflag = 1;
%                     else
%                         if currsort == length(uvsorts)
%                             snpidx2 = size(tmpv2,2);
%                         else
%                             snpidx2 = clusize(currsort);
%                         end
%                         partialflag = 0;
%                     end
%                 end


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


% --- Executes on button press in skip.
function skip_Callback(hObject, eventdata, handles)
% hObject    handle to skip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv tmpv2 goodflag r5pc s5pc
global heightsmean heightsstd widthsmean widthsstd ens tmpclusters fname retry tmpsortidx tmpunsortedidx vs tridx tmps2 clusize undohistory

bdf2 = disable_all(handles);

if max(1, length(unique(tmpsortidx(tmpsortidx>0))))==1
    
    PF2 = PF;
    tmpsortidx(tmpunsortedidx) = 1;

    for vt = validtrials
            PF2.rec(vt).v_final{1} = tmpv2(:,(tridx== vt));
            PF2.rec(vt).plexon_spikes_final{1} = tmps2(tridx== vt);
    end


    [rasters, FrameRate_Hz] = tmpraster(PF2);

    r1 = compress_raster(rasters{1});

    sc=calcFstatistics2(r1);
    
    for vt = validtrials

        PF.rec(vt).v_MU{1,ch} = tmpv2(:,((tridx== vt)));
        PF.rec(vt).plexon_spikes_MU{1,ch} = tmps2(((tridx== vt)));
        PF.rec(vt).F_MU{1,ch} = sc;

    end

    ens{ch}=1;

else
    if goodflag==1
        fc = finalcluster;
        uiwait(fc)
    end
end
goodflag = 0;

if ~retry
    ch = ch + 1;
end

if ch >= 97
%     load(fname)
     s1 = strfind(fname, '/');
% 
%     for vt = validtrials
%         PF.rec(vt).v_final = tmpclusters(vt).v;
%         PF.rec(vt).plexon_spikes_final = tmpclusters(vt).plexon_spikes;
%     end

    save([fname(1:end-3) 'rs.mat'], 'PF')

    save([fname(1:s1(end)) 'sortparams.mat'], 'ens', 'heightsmean', 'heightsstd', 'widthsmean', 'widthsstd', 'r5pc', 's5pc')
    fprintf('Done')
    return
end

[vsorts, vtrials] = find(cellfun(@(x) ~isempty(x),vv(:,ch,:))==1);

uvsorts = unique(vsorts);

% uvidx = 0;
% clusize = [];
% 
% for uv = uvsorts'
%     tmpv = cat(2,vv{uv, ch, :});
%     tmpv2(:,uvidx+1:uvidx+size(tmpv,2)) = tmpv;
%     uvidx = uvidx+size(tmpv,2);
%     clusize = [clusize uvidx];
% end

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

set(handles.ChanNum, 'String', ['Channel: ' num2str(ch)]);


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
undohistory{1}.currsort = currsort;

% tmpsortidx = zeros(1, size(tmpv2,2));
% tmpunsortedidx = find(tmpsortidx==0);

snpidx = intersect(snpidx1:snpidx2, tmpunsortedidx);

% tmpv2 = cat(2,vv{uvsorts, ch, :});
    
%     if mean(tmpv2(1,:))<mean(tmpv2(7,:))
% 
%         continue
%     else
        goodflag=1;
%     end


    [plotmin, plotmax] = findplotminmax(tmpv2);
    allclusters

%     tmpv2 = cat(2,vv{uvsorts, ch, :});

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



plot(handles.mainplot, tmpv2(:,snpidx)); set(handles.mainplot, 'YLim', [plotmin plotmax]);

kms=6;
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


% --- Executes on button press in undo.
function undo_Callback(hObject, eventdata, handles)
% hObject    handle to undo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels validtrials tmpv tmpv2 tmpclusters r5pc s5pc 
global heightsmean heightsstd widthsmean widthsstd goodflag fname tmpsortidx tmpunsortedidx clusize vs tridx tmps2 retry undohistory

% tmpsortidx(snpidx(find(labels==labelnum))) = clu;

if length(undohistory)>1
    tmpsortidx = undohistory{end-1}.tmpsortidx;
    snpidx1 = undohistory{end}.snpidx1;
    snpidx2 = undohistory{end}.snpidx2;
    currsort = undohistory{end}.currsort;
    partialflag = undohistory{end}.partialflag;
    undohistory(end) = [];

    tmpunsortedidx = find(tmpsortidx==0);
end
length(undohistory)
tempcluster
set(handles.figure1, 'Visible', 'on')

kms = 6;

set(handles.ChanNum, 'String', ['Channel: ' num2str(ch)]);


snpidx = intersect(snpidx1:snpidx2, tmpunsortedidx);

if ~isempty(snpidx)
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
end


% --- Executes on button press in nextclusterbutton.
function nextclusterbutton_Callback(hObject, eventdata, handles)
% hObject    handle to nextclusterbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in skipSU.
function skipSU_Callback(hObject, eventdata, handles)
% hObject    handle to skipSU (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over skipSU.
function skipSU_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to skipSU (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
