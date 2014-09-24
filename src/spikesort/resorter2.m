function varargout = resorter2(varargin)
% RESORTER2 MATLAB code for resorter2.fig
%      RESORTER2, by itself, creates a new RESORTER2 or raises the existing
%      singleton*.
%
%      H = RESORTER2 returns the handle to a new RESORTER2 or the handle to
%      the existing singleton*.
%
%      RESORTER2('CALLBACK',hObject,eventData,handles,...) calls the localff
%      function named CALLBACK in RESORTER2.M with the given input arguments.
%
%      RESORTER2('Property','Value',...) creates a new RESORTER2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before resorter2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to resorter2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help resorter2

% Last Modified by GUIDE v2.5 27-Mar-2012 18:32:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @resorter2_OpeningFcn, ...
                   'gui_OutputFcn',  @resorter2_OutputFcn, ...
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


% --- Executes just before resorter2 is made visible.
function resorter2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to resorter2 (see VARARGIN)

% Choose default command line output for resorter2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes resorter2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = resorter2_OutputFcn(hObject, eventdata, handles) 
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
set(handles.status, 'String', 'Wait'); drawnow
set(handles.nextbutton, 'Enable', 'off')

addpath /auto/k1/moliver/code/kdica
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

global partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv goodflag tmpclusters fname sortparams
global heightsmean heightsstd widthsmean widthsstd s5pc r5pc ens
s1 = strfind(fname, '/');

if exist([fname(1:s1(end)) 'sortparams.mat'])
    sortparams = load([fname(1:s1(end)) 'sortparams.mat']);
    ens = sortparams.ens;
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

        csort = uvsorts(currsort);

        tmpv = cat(2,vv{csort, ch, vtrials(vsorts==csort)});
        tmpv2 = cat(2,vv{uvsorts, ch, :});
        plotmin = min(tmpv2(:));
        plotmax = max(tmpv2(:));
        
        plotokflag = 0;
        
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
        end
        
        tmpv2 = cat(2,vv{uvsorts, ch, :});
            
        
        
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

        if size(tmpv,2) > str2num(get(handles.snipnum, 'string'));
            snpidx1 = 1;
            snpidx2 = str2num(get(handles.snipnum, 'string'));
            partialflag = 1;
        else
            snpidx1 = 1;
            snpidx2 = size(tmpv,2);
            partialflag = 0;
        end

        snpidx = snpidx1:snpidx2;

        plot(handles.mainplot, tmpv(:,snpidx)); set(handles.mainplot, 'YLim', [plotmin plotmax]);


        if get(handles.slope, 'Value')==1
            labels = cluster(diff(tmpv(:,snpidx))', kms, handles);
        else
            labels = cluster(tmpv(:,snpidx)', kms, handles);
        end



        for kk = 1:kms
            hh = ['handles.kmean' num2str(kk)];
            bdf = get(eval(hh), 'ButtonDownFcn');
            cla(eval(hh), 'reset')
            plot(eval(hh), tmpv(:,snpidx(labels==kk))); set(eval(hh), 'YLim', [plotmin plotmax]);
            set(eval(hh), 'ButtonDownFcn', bdf)
        end

        flag = 1;
    end
end
set(handles.status, 'String', 'Ready');
set(handles.nextbutton, 'Enable', 'on');

% --- Executes on button press in nextbutton.
function nextbutton_Callback(hObject, eventdata, handles)
% hObject    handle to nextbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv tmpv2 goodflag r5pc s5pc
global heightsmean heightsstd widthsmean widthsstd ens tmpclusters fname
set(handles.status, 'String', 'Wait'); drawnow
set(handles.nextbutton, 'Enable', 'off')
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
            ch = ch + 1;
            
            if ch >= 97
                continue
            end
            
            currsort = 1;
            [vsorts, vtrials] = find(cellfun(@(x) ~isempty(x),vv(:,ch,:))==1);

            uvsorts = unique(vsorts);

        end

        set(handles.ChanNum, 'String', ['Channel: ' num2str(ch)]);

        csort = uvsorts(currsort);

        tmpv = cat(2,vv{csort, ch, vtrials(vsorts==csort)});

        
        if size(tmpv,2) > str2num(get(handles.snipnum, 'string'));
            snpidx1 = 1;
            snpidx2 = str2num(get(handles.snipnum, 'string'));
            partialflag = 1;
        else
            snpidx1 = 1;
            snpidx2 = size(tmpv,2);
            partialflag = 0;
        end
        
        
    elseif partialflag == 1 && ch < 97
        
        snpidx1 = snpidx1 + str2num(get(handles.snipnum, 'string'));
        snpidx2 = snpidx2 + str2num(get(handles.snipnum, 'string'));
        
        if size(tmpv,2) < snpidx2
            snpidx2 = size(tmpv,2);
            partialflag = 0;
        else
            partialflag = 1;
        end
        
    else
        load(fname)
        s1 = strfind(fname, '/');
        
        for vt = validtrials
            PF.rec(vt).v_final = tmpclusters(vt).v;
            PF.rec(vt).plexon_spikes_final = tmpclusters(vt).plexon_spikes;
        end
        
        save(fname, 'PF')
        
        save([fname(1:s1(end)) 'sortparams.mat'], 'ens', 'heightsmean', 'heightsstd', 'widthsmean', 'widthsstd', 'r5pc', 's5pc','tmpclusters')
        return
    end
            
    
    snpidx = snpidx1:snpidx2;
    tmpv2 = cat(2,vv{uvsorts, ch, :});
    
%     if mean(tmpv2(1,:))<mean(tmpv2(7,:))
% 
%         continue
%     else
        goodflag=1;
%     end

    if (currsort == 1 && snpidx1 == 1) || isempty(heightsmean)
%         keyboard
        plotmin = min(tmpv2(:));
        plotmax = max(tmpv2(:));
        
        plotokflag = 0;
        
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
        end
        
        tmpv2 = cat(2,vv{uvsorts, ch, :});
        

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

    
    plot(handles.mainplot, tmpv(:,snpidx)); set(handles.mainplot, 'YLim', [plotmin plotmax]);

    
	if get(handles.slope, 'Value')==1
        labels = cluster(diff(tmpv(:,snpidx))', kms, handles);
    else
        labels = cluster(tmpv(:,snpidx)', kms, handles);
    end

    for kk = 1:kms
        hh = ['handles.kmean' num2str(kk)];
        bdf = get(eval(hh), 'ButtonDownFcn');
        cla(eval(hh), 'reset')
        plot(eval(hh), tmpv(:,snpidx(labels==kk))); set(eval(hh), 'YLim', [plotmin plotmax]);
        set(eval(hh), 'ButtonDownFcn', bdf)
    end
        
    flag = 1;
end

set(handles.status, 'String', 'Ready');
set(handles.nextbutton, 'Enable', 'on');
% --- Executes on button press in deletesc.
function deletesc_Callback(hObject, eventdata, handles)
% hObject    handle to deletesc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of deletesc
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
global partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv2
set(handles.status, 'String', 'Wait'); drawnow
set(handles.nextbutton, 'Enable', 'off')
vv = cat(3,PF.rec(:).v);

kms = 6;

set(handles.ChanNum, 'String', ['Channel: ' num2str(ch)]);

csort = uvsorts(currsort);

tmpv = cat(2,vv{csort, ch, vtrials(vsorts==csort)});
tmpv2 = cat(2,vv{uvsorts, ch, :});
% keyboard
plot(handles.mainplot, tmpv(:,snpidx)); set(handles.mainplot, 'YLim', [plotmin plotmax]);
% keyboard

if get(handles.slope, 'Value')==1
    labels = cluster(diff(tmpv(:,snpidx))', kms, handles);
else
    labels = cluster(tmpv(:,snpidx)', kms, handles);
end


for kk = 1:kms
    hh = ['handles.kmean' num2str(kk)];
    bdf = get(eval(hh), 'ButtonDownFcn');
    cla(eval(hh), 'reset')
    plot(eval(hh), tmpv(:,snpidx(labels==kk))); set(eval(hh), 'YLim', [plotmin plotmax]);
    set(eval(hh), 'ButtonDownFcn', bdf)
end
set(handles.status, 'String', 'Ready');
set(handles.nextbutton, 'Enable', 'on');


% --- Executes on button press in peaks.
function pca_Callback(hObject, eventdata, handles)
% hObject    handle to peaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global r5pc s5pc partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv2
set(handles.status, 'String', 'Wait'); drawnow
set(handles.nextbutton, 'Enable', 'off')
vv = cat(3,PF.rec(:).v);

kms = 6;

set(handles.ChanNum, 'String', ['Channel: ' num2str(ch)]);

csort = uvsorts(currsort);

tmpv = cat(2,vv{csort, ch, vtrials(vsorts==csort)});
tmpv2 = cat(2,vv{uvsorts, ch, :});

plot(handles.mainplot, tmpv(:,snpidx)); set(handles.mainplot, 'YLim', [plotmin plotmax]);

if get(handles.slope, 'Value')==1

    if str2num(get(handles.pcnum, 'String'))>0
        ktmpv = (s5pc{ch}(:,str2num(get(handles.pcnum, 'String')))'*diff(tmpv(:,snpidx)));
    else
        ktmpv = (s5pc{ch}(:,1)'*diff(tmpv(:,snpidx)));
    end
else
    if str2num(get(handles.pcnum, 'String'))>0
        ktmpv = (r5pc{ch}(:,str2num(get(handles.pcnum, 'String')))'*tmpv(:,snpidx));
    else
        ktmpv = (s5pc{ch}(:,1)'*diff(tmpv(:,snpidx)));
    end
end

labels = cluster(ktmpv', kms, handles);


for kk = 1:kms
    hh = ['handles.kmean' num2str(kk)];
    bdf = get(eval(hh), 'ButtonDownFcn');
    cla(eval(hh), 'reset')
    plot(eval(hh), tmpv(:,snpidx(labels==kk))); set(eval(hh), 'YLim', [plotmin plotmax]);
    set(eval(hh), 'ButtonDownFcn', bdf)
end
set(handles.status, 'String', 'Ready');
set(handles.nextbutton, 'Enable', 'on');


% --- Executes on button press in peaks.
function peaks_Callback(hObject, eventdata, handles)
% hObject    handle to peaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global r5pc s5pc partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv2
global heightsmean heightsstd widthsmean widthsstd
set(handles.status, 'String', 'Wait'); drawnow
set(handles.nextbutton, 'Enable', 'off')
vv = cat(3,PF.rec(:).v);

kms = 6;

set(handles.ChanNum, 'String', ['Channel: ' num2str(ch)]);

csort = uvsorts(currsort);

tmpv = cat(2,vv{csort, ch, vtrials(vsorts==csort)});
tmpv2 = cat(2,vv{uvsorts, ch, :});

plot(handles.mainplot, tmpv(:,snpidx)); set(handles.mainplot, 'YLim', [plotmin plotmax]);

heights = max(tmpv(:,snpidx)) - min(tmpv(:,snpidx));
heights = heights - heightsmean;
heights = heights ./ heightsstd;
[minr,minc] = find(bsxfun(@(x,y) (x == y), tmpv(:,snpidx), min(tmpv(:,snpidx)))==1);
[maxr,maxc] = find(bsxfun(@(x,y) (x == y), tmpv(:,snpidx), max(tmpv(:,snpidx)))==1);
minr = minr(unique(minc));
maxr = maxr(unique(maxc));

widths = maxr - minr;
widths = widths - widthsmean;
widths = widths ./ widthsstd;

% keyboard

if get(handles.slope, 'Value')==1
    if str2num(get(handles.pcnum, 'String'))>0
        ktmpv2 = (s5pc{ch}(:,str2num(get(handles.pcnum, 'String')))'*diff(tmpv(:,snpidx)));
    else
        ktmpv2 = [];
    end
else
    if str2num(get(handles.pcnum, 'String'))>0
        ktmpv2 = (r5pc{ch}(:,str2num(get(handles.pcnum, 'String')))'*tmpv(:,snpidx));
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
    plot(eval(hh), tmpv(:,snpidx(labels==kk))); set(eval(hh), 'YLim', [plotmin plotmax]);
    set(eval(hh), 'ButtonDownFcn', bdf)
end
set(handles.status, 'String', 'Ready');
set(handles.nextbutton, 'Enable', 'on');


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


function labels = cluster(data, k, handles)

if get(handles.spectralcluster, 'Value') == 1

    W = SimGraph(data, 1);
    A = SpectralClustering(W, k, 3);
    labels = full(sum(bsxfun(@times, A, 1:6),2));
elseif get(handles.kmeanscluster, 'Value') == 1
    
    if size(data,1)<size(data,2)
        k=min(k, size(data,1));
        labels = fkmeans(data,k);
    else
        labels = kmeans(data,k);
    end
elseif get(handles.kmedoidscluster, 'Value') == 1
    
    if size(data,1)==1
        labels=1;
    else
        k=min(k, size(data,1));
        labels = kmedoids(data', k)';
    end
end


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
global partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv tmpv2 tmpclusters r5pc s5pc 
global heightsmean heightsstd widthsmean widthsstd goodflag fname ens
set(handles.status, 'String', 'Wait');  drawnow
set(handles.nextbutton, 'Enable', 'off')
if get(handles.deletesc, 'Value') == 1

    [tvidx] = snpidx(find(labels==labelnum));
    
    trsz = cat(1,squeeze(cellfun(@(x) size(x,2), vv(csort, ch, vtrials(vsorts==csort)))));
    cnt = 0;
    for vt = vtrials(vsorts==csort)'
        cnt = cnt+1;
        delidx = tvidx(tvidx <= trsz(cnt));
        tvidx = tvidx(tvidx > trsz(cnt));
        tvidx = tvidx - trsz(cnt);
        PF.rec(validtrials(vt)).v{csort,ch}(:,delidx) = [];
        PF.rec(validtrials(vt)).plexon_spikes{csort,ch}(delidx) = []; 
    end
    clu = 0;
elseif get(handles.cluster1, 'Value') == 1
    
    clu=1;
elseif get(handles.cluster2, 'Value') == 1
    
    clu=2;
elseif get(handles.cluster3, 'Value') == 1
    
    clu=3;
elseif get(handles.cluster4, 'Value') == 1
    
    clu=4;
elseif get(handles.cluster5, 'Value') == 1
    
    clu=5;
elseif get(handles.cluster6, 'Value') == 1
    
    clu=6;
else
    clu=0;
end

if clu>0
    [tvidx] = snpidx(find(labels==labelnum));
    
%     trsz = cat(1,squeeze(cellfun(@(x) size(x,2), vv(csort, ch, vtrials(vsorts==csort)))));

    trsz = cat(1,squeeze(cellfun(@(x) size(x,2), vv(csort, ch, vtrials))));
    trsz(vsorts~=csort)=0;
    cnt = 0;
    
%     for vt = vtrials(vsorts==csort)'
    for vt = vtrials'
        cnt = cnt+1;
        delidx = tvidx(tvidx <= trsz(cnt));
        tvidx = tvidx(tvidx > trsz(cnt));
        tvidx = tvidx - trsz(cnt);
        if size(tmpclusters(validtrials(vt)).v,1) < clu
            tmpclusters(validtrials(vt)).v{clu,ch} = [];
            tmpclusters(validtrials(vt)).plexon_spikes{clu,ch} = [];
        end
        if size(tmpclusters(validtrials(vt)).v,2) < ch
            tmpclusters(validtrials(vt)).v{clu,ch} = [];
            tmpclusters(validtrials(vt)).plexon_spikes{clu,ch} = [];
        end
        
        tmpclusters(validtrials(vt)).v{clu,ch} = [tmpclusters(validtrials(vt)).v{clu,ch} PF.rec(validtrials(vt)).v{csort,ch}(:,delidx)];
        tmpclusters(validtrials(vt)).plexon_spikes{clu,ch} = [tmpclusters(validtrials(vt)).plexon_spikes{clu,ch} PF.rec(validtrials(vt)).plexon_spikes{csort,ch}(delidx)];
        [~,dsidx] = sort(tmpclusters(validtrials(vt)).plexon_spikes{clu,ch}, 'ascend');
        tmpclusters(validtrials(vt)).v{clu,ch} = tmpclusters(validtrials(vt)).v{clu,ch}(:,dsidx);
        tmpclusters(validtrials(vt)).plexon_spikes{clu,ch} = tmpclusters(validtrials(vt)).plexon_spikes{clu,ch}(dsidx);
        PF.rec(validtrials(vt)).v{csort,ch}(:,delidx) = [];
        PF.rec(validtrials(vt)).plexon_spikes{csort,ch}(delidx) = [];
    end
    tempcluster
    set(handles.figure1, 'Visible', 'on')
end

vv = cat(3,PF.rec(:).v);

kms = 6;

set(handles.ChanNum, 'String', ['Channel: ' num2str(ch)]);

csort = uvsorts(currsort);

tmpv = cat(2,vv{csort, ch, vtrials(vsorts==csort)});
tmpv2 = cat(2,vv{uvsorts, ch, :});


if size(tmpv,2) < snpidx2;
    snpidx2 = size(tmpv,2);
    snpidx = snpidx1:snpidx2;
    partialflag = 0;
end

if ~isempty(tmpv)
%     if isempty(tmpv(:,snpidx))
%         keyboard
%     end
    
    plot(handles.mainplot, tmpv(:,snpidx)); set(handles.mainplot, 'YLim', [plotmin plotmax]);

    if get(handles.slope, 'Value')==1
        labels = cluster(diff(tmpv(:,snpidx))', kms, handles);
    else
        labels = cluster(tmpv(:,snpidx)', kms, handles);
    end


    for kk = 1:kms
        hh = ['handles.kmean' num2str(kk)];
        bdf = get(eval(hh), 'ButtonDownFcn');
        cla(eval(hh), 'reset')
        plot(eval(hh), tmpv(:,snpidx(labels==kk))); set(eval(hh), 'YLim', [plotmin plotmax]);
        set(eval(hh), 'ButtonDownFcn', bdf)
    end
else
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
                ch = ch + 1;
                
                if ch >= 97
                    continue
                end

                currsort = 1;
                [vsorts, vtrials] = find(cellfun(@(x) ~isempty(x),vv(:,ch,:))==1);

                uvsorts = unique(vsorts);

            end

            set(handles.ChanNum, 'String', ['Channel: ' num2str(ch)]);

            csort = uvsorts(currsort);

            tmpv = cat(2,vv{csort, ch, vtrials(vsorts==csort)});


            if size(tmpv,2) > str2num(get(handles.snipnum, 'string'));
                snpidx1 = 1;
                snpidx2 = str2num(get(handles.snipnum, 'string'));
                partialflag = 1;
            else
                snpidx1 = 1;
                snpidx2 = size(tmpv,2);
                partialflag = 0;
            end


        elseif partialflag == 1 && ch < 97

            snpidx1 = snpidx1 + str2num(get(handles.snipnum, 'string'));
            snpidx2 = snpidx2 + str2num(get(handles.snipnum, 'string'));

            if size(tmpv,2) < snpidx2
                snpidx2 = size(tmpv,2);
                partialflag = 0;
            else
                partialflag = 1;
            end
        else
            load(fname)
            s1 = strfind(fname, '/');

            for vt = validtrials
                PF.rec(vt).v_final = tmpclusters(vt).v;
                PF.rec(vt).plexon_spikes_final = tmpclusters(vt).plexon_spikes;
            end

            save(fname, 'PF')

            save([fname(1:s1(end)) 'sortparams.mat'], 'ens', 'heightsmean', 'heightsstd', 'widthsmean', 'widthsstd', 'r5pc', 's5pc', 'tmpclusters')
            fprintf('Done')
            return
            
        end


        snpidx = snpidx1:snpidx2;
        tmpv2 = cat(2,vv{uvsorts, ch, :});

%         if mean(tmpv2(1,:))<mean(tmpv2(7,:))
% 
%             continue
%         else
            goodflag=1;
%         end

        if (currsort == 1 && snpidx1 == 1) || isempty(heightsmean)
            plotmin = min(tmpv2(:));
            plotmax = max(tmpv2(:));

            plotokflag = 0;
        
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
            end
            
            tmpv2 = cat(2,vv{uvsorts, ch, :});
            
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


        plot(handles.mainplot, tmpv(:,snpidx)); set(handles.mainplot, 'YLim', [plotmin plotmax]);

        if get(handles.slope, 'Value')==1
            labels = cluster(diff(tmpv(:,snpidx))', kms, handles);
        else
            labels = cluster(tmpv(:,snpidx)', kms, handles);
        end

        for kk = 1:kms
            hh = ['handles.kmean' num2str(kk)];
            bdf = get(eval(hh), 'ButtonDownFcn');
            cla(eval(hh), 'reset')
            plot(eval(hh), tmpv(:,snpidx(labels==kk))); set(eval(hh), 'YLim', [plotmin plotmax]);
            set(eval(hh), 'ButtonDownFcn', bdf)
        end

        flag = 1;
    end
end
set(handles.status, 'String', 'Ready');
set(handles.nextbutton, 'Enable', 'on');


% --- Executes on button press in skip.
function skip_Callback(hObject, eventdata, handles)
% hObject    handle to skip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv tmpv2 goodflag r5pc s5pc
global heightsmean heightsstd widthsmean widthsstd ens tmpclusters fname

set(handles.status, 'String', 'Wait');  drawnow
set(handles.nextbutton, 'Enable', 'off')
set(handles.skip, 'Enable', 'off')

if goodflag==1
    fc = finalcluster;
    uiwait(fc)
end
goodflag = 0;
ch = ch + 1;

if ch >= 97
    load(fname)
    s1 = strfind(fname, '/');

    for vt = validtrials
        PF.rec(vt).v_final = tmpclusters(vt).v;
        PF.rec(vt).plexon_spikes_final = tmpclusters(vt).plexon_spikes;
    end

    save(fname, 'PF')

    save([fname(1:s1(end)) 'sortparams.mat'], 'ens', 'heightsmean', 'heightsstd', 'widthsmean', 'widthsstd', 'r5pc', 's5pc', 'tmpclusters')
    fprintf('Done')
    return
end

currsort = 1;
[vsorts, vtrials] = find(cellfun(@(x) ~isempty(x),vv(:,ch,:))==1);

uvsorts = unique(vsorts);

set(handles.ChanNum, 'String', ['Channel: ' num2str(ch)]);

csort = uvsorts(currsort);

tmpv = cat(2,vv{csort, ch, vtrials(vsorts==csort)});


if size(tmpv,2) > str2num(get(handles.snipnum, 'string'));
    snpidx1 = 1;
    snpidx2 = str2num(get(handles.snipnum, 'string'));
    partialflag = 1;
else
    snpidx1 = 1;
    snpidx2 = size(tmpv,2);
    partialflag = 0;
end


snpidx = snpidx1:snpidx2;
tmpv2 = cat(2,vv{uvsorts, ch, :});
    
%     if mean(tmpv2(1,:))<mean(tmpv2(7,:))
% 
%         continue
%     else
        goodflag=1;
%     end


    plotmin = min(tmpv2(:));
    plotmax = max(tmpv2(:));

    plotokflag = 0;

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
    end

    tmpv2 = cat(2,vv{uvsorts, ch, :});

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



plot(handles.mainplot, tmpv(:,snpidx)); set(handles.mainplot, 'YLim', [plotmin plotmax]);

kms=6;
if get(handles.slope, 'Value')==1
    labels = cluster(diff(tmpv(:,snpidx))', kms, handles);
else
    labels = cluster(tmpv(:,snpidx)', kms, handles);
end

for kk = 1:kms
    hh = ['handles.kmean' num2str(kk)];
    bdf = get(eval(hh), 'ButtonDownFcn');
    cla(eval(hh), 'reset')
    plot(eval(hh), tmpv(:,snpidx(labels==kk))); set(eval(hh), 'YLim', [plotmin plotmax]);
    set(eval(hh), 'ButtonDownFcn', bdf)
end

flag = 1;


set(handles.status, 'String', 'Ready');
set(handles.nextbutton, 'Enable', 'on');
set(handles.skip, 'Enable', 'on');
