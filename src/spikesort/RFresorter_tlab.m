function varargout = RFresorter_tlab(varargin)
% RFRESORTER_TLAB MATLAB code for RFresorter_tlab.fig
%      RFRESORTER_TLAB, by itself, creates a new RFRESORTER_TLAB or raises the existing
%      singleton*.
%
%      H = RFRESORTER_TLAB returns the handle to a new RFRESORTER_TLAB or the handle to
%      the existing singleton*.
%
%      RFRESORTER_TLAB('CALLBACK',hObject,eventData,handles,...) calls the localff
%      function named CALLBACK in RFRESORTER_TLAB.M with the given input arguments.
%
%      RFRESORTER_TLAB('Property','Value',...) creates a new RFRESORTER_TLAB or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RFresorter_tlab_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RFresorter_tlab_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RFresorter_tlab

% Last Modified by GUIDE v2.5 13-Mar-2013 23:27:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RFresorter_tlab_OpeningFcn, ...
                   'gui_OutputFcn',  @RFresorter_tlab_OutputFcn, ...
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


% --- Executes just before RFresorter_tlab is made visible.
function RFresorter_tlab_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RFresorter_tlab (see VARARGIN)

% Choose default command line output for RFresorter_tlab
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% sort_function(handles)
% UIWAIT makes RFresorter_tlab wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RFresorter_tlab_OutputFcn(hObject, eventdata, handles) 
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


% --- Executes after reading file or on retry.
function sort_function(handles)
% hObject    handle to sortbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
bdf2 = disable_all(handles);
% global uvsorts
% clear uvsorts
% clear PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv tmpv2 goodflag fname
global fname 
set(handles.filename, 'String', fname);

global snips partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv goodflag tmpclusters sortparams
global heightsmean heightsstd widthsmean widthsstd s5pc r5pc ens tmpsortidx clusize tmpv2 tmpunsortedidx vs tridx tmps2 retry undohistory colortable resorter_handles
global msnips mdsnips

uvsorts = [];
resorter_handles = handles;


retry=0;


        kms = 6;


       
        tmpsortidx = zeros(1, size(snips,1));
        tmpunsortedidx = find(tmpsortidx==0);
  
        % PCA for slopes
        mdsnips = mean(diff(snips')');
        s5pcsnips = bsxfun(@minus, diff(snips')', mdsnips);
        covtv = s5pcsnips'*s5pcsnips;
        [U,S,V]= svd(covtv);

        s5pc = U;

        % PCA for snips
        msnips = mean(snips);
        r5pcsnips = bsxfun(@minus, snips, msnips);
        covtv = r5pcsnips'*r5pcsnips;
        [U,S,V]= svd(covtv);

        r5pc = U;

        heights = max(snips') - min(snips');
        heightsmean = mean(heights);
        heights = heights - heightsmean;
        heightsstd = std(heights);
        [minr,minc] = find(bsxfun(@(x,y) (x == y), snips', min(snips'))==1);
        [maxr,maxc] = find(bsxfun(@(x,y) (x == y), snips', max(snips'))==1);
        minr = minr(unique(minc));
        maxr = maxr(unique(maxc));
        widths = minr - maxr;
        widthsmean = mean(widths);
        widths = widths - widthsmean;
        widthsstd = std(widths);

        
        [plotmin, plotmax] = findplotminmax(snips);

        if size(snips,1) > str2num(get(handles.snipnum, 'string'));
            snpidx1 = 1;
            snpidx2 = str2num(get(handles.snipnum, 'string'));
            partialflag = 1;
        else
            snpidx1 = 1;
            snpidx2 = size(snips,1);
            partialflag = 0;
        end
        
        undohistory{1}.tmpsortidx=tmpsortidx;
        undohistory{1}.snpidx1=snpidx1;
        undohistory{1}.snpidx2=snpidx2;
        undohistory{1}.currsort=currsort;
        undohistory{1}.partialflag = partialflag;
        snpidx = snpidx1:snpidx2;


        if get(handles.slope, 'Value')==1
            labels = cluster(diff(snips(snpidx,:)')', kms, handles);
        else
            labels = cluster(snips(snpidx,:), kms, handles);
        end

        plot_main(handles);



enable_all(handles, bdf2)

% --- Executes on button press in nextbutton.
function nextbutton_Callback(hObject, eventdata, handles)
% hObject    handle to nextbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv tmpv2 goodflag r5pc s5pc shuffle
global heightsmean snips heightsstd widthsmean widthsstd ens tmpclusters fname tmpsortidx clusize tmpunsortedidx vs tridx tmps2 retry undohistory skip_button
bdf2 = disable_all(handles);
% if isempty(ch)
%     ch = 1;
% end

if isempty(currsort)
    currsort=0;
end

flag = 0;
kms = 6;

     
%     snpidx1 = snpidx1 + str2num(get(handles.snipnum, 'string'));
%     snpidx2 = snpidx2 + str2num(get(handles.snipnum, 'string'));

     snpidx1 = snpidx2+1;
     snpidx2 = snpidx1 + str2num(get(handles.snipnum, 'string')) -1;

    if size(snips,1) <= snpidx2
        snpidx2 = size(snips,1);
        partialflag = 0;
    else
        partialflag = 1;
    end
        
   
    
    snpidx = intersect(shuffle(snpidx1:snpidx2), tmpunsortedidx);
 
    

    if get(handles.slope, 'Value')==1
        labels = cluster(diff(snips(snpidx,:)')', kms, handles);
    else
        labels = cluster(snips(snpidx,:), kms, handles);
    end
    
    
    plot_main(handles);

        
    flag = 1;
% end

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
global partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv2 vv snips
bdf2 = disable_all(handles);
% vv = cat(3,PF.rec(:).v);

kms = 6;


if get(handles.slope, 'Value')==1
    labels = cluster(diff(snips(snpidx,:)')', kms, handles);
else
    labels = cluster(snips(snpidx,:), kms, handles);
end

plot_main(handles);

enable_all(handles, bdf2)


% --- Executes on button press in peaks.
function pca_Callback(hObject, eventdata, handles)
% hObject    handle to peaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global r5pc s5pc partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv2 vv snips
global msnips mdsnips

bdf2 = disable_all(handles);

kms = 6;

if get(handles.slope, 'Value')==1
    dsnips = bsxfun(@minus, diff(snips(snpidx,:)')', mdsnips);
    if str2num(get(handles.pcnum, 'String'))>0
        ktmpv = (s5pc(:,1:str2num(get(handles.pcnum, 'String')))'*dsnips');
    else
        ktmpv = (s5pc(:,1:5)'*diff(snips(snpidx,:)'));
    end
else
    dsnips = bsxfun(@minus, snips(snpidx,:), msnips);
    if str2num(get(handles.pcnum, 'String'))>0
        ktmpv = (r5pc(:,1:str2num(get(handles.pcnum, 'String')))'*dsnips');
    else
        ktmpv = (r5pc(:,1:5)'*dsnips');
    end
end

labels = cluster(ktmpv', kms, handles);

plot_main(handles);

enable_all(handles, bdf2)


% --- Executes on button press in peaks.
function peaks_Callback(hObject, eventdata, handles)
% hObject    handle to peaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global r5pc s5pc partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv2 snips
global heightsmean heightsstd widthsmean widthsstd
global msnips mdsnips

bdf2 = disable_all(handles);


kms = 6;

heights = max(snips(snpidx,:)') - min(snips(snpidx,:)');
heights = heights - heightsmean;
heights = heights ./ heightsstd;
[minr,minc] = find(bsxfun(@(x,y) (x == y), snips(snpidx,:)', min(snips(snpidx,:)'))==1);
[maxr,maxc] = find(bsxfun(@(x,y) (x == y), snips(snpidx,:)', max(snips(snpidx,:)'))==1);
minr = minr(unique(minc));
maxr = maxr(unique(maxc));

widths = maxr - minr;
widths = widths - widthsmean;
widths = widths ./ widthsstd;


if get(handles.slope, 'Value')==1
    dsnips = bsxfun(@minus, diff(snips(snpidx,:)')', mdsnips);
    if str2num(get(handles.pcnum, 'String'))>0
        ktmpv2 = (s5pc(:,1:str2num(get(handles.pcnum, 'String')))'*dsnips');
    else
        ktmpv2 = [];
    end
else
    dsnips = bsxfun(@minus, snips(snpidx,:), msnips);
    if str2num(get(handles.pcnum, 'String'))>0
        ktmpv2 = (r5pc(:,1:str2num(get(handles.pcnum, 'String')))'*dsnips');
    else
        ktmpv2 = [];
    end
end


ktmpv = [widths'; heights; ktmpv2];


labels = cluster(ktmpv', kms, handles);

plot_main(handles);

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
global partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels validtrials tmpv tmpv2 r5pc s5pc snips shuffle
global heightsmean heightsstd widthsmean widthsstd goodflag fname tmpsortidx tmpunsortedidx clusize vs tridx tmps2 retry undohistory clu ens

bdf2 = disable_all(handles);

if isempty(clu); clu=1;end
tmpsortidx(snpidx(find(labels==labelnum))) = clu;
tmpunsortedidx = find(tmpsortidx==0);

if clu>0

    tempcluster_tlab;
    set(handles.figure1, 'Visible', 'on');
end

% vv = cat(3,PF.rec(:).v);

kms = 6;

% set(handles.ChanNum, 'String', ['Channel: ' num2str(ch)]);

undohistory{end+1}.tmpsortidx = tmpsortidx;
undohistory{end}.snpidx1 = snpidx1;
undohistory{end}.snpidx2 = snpidx2;
undohistory{end}.currsort=currsort;
undohistory{end}.partialflag = partialflag;


snpidx = intersect(shuffle(snpidx1:snpidx2), tmpunsortedidx);
    


if ~isempty(snips(snpidx,:)')

    


    if get(handles.slope, 'Value')==1
        labels = cluster(diff(snips(snpidx,:)')', kms, handles);
    else
        labels = cluster(snips(snpidx,:), kms, handles);
    end
    
    plot_main(handles);

else


end
enable_all(handles, bdf2)




function bdf = disable_all(handles)
set(handles.status, 'String', 'Wait');  drawnow
set(handles.nextbutton, 'Enable', 'off')
% set(handles.skip, 'Enable', 'off')
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
% set(handles.skip, 'Enable', 'on');
set(handles.kmean1, 'ButtonDownFcn', bdf{1});
set(handles.kmean2, 'ButtonDownFcn', bdf{2});
set(handles.kmean3, 'ButtonDownFcn', bdf{3});
set(handles.kmean4, 'ButtonDownFcn', bdf{4});
set(handles.kmean5, 'ButtonDownFcn', bdf{5});
set(handles.kmean6, 'ButtonDownFcn', bdf{6});

function plot_main(handles)
global snpidx plotmin plotmax snips r5pc
global msnips mdsnips labels

kms = 6;   % the number of clusters
colorvals = ['y', 'm', 'c', 'r', 'g', 'b'];

% First draw the main graph
if (get(handles.Traces, 'Value') == 1)
    for kk=1:kms
        plot(handles.mainplot, snips(snpidx(labels==kk),:)','Color', colorvals(kk));
        hold(handles.mainplot,'on');
    end
    set(handles.mainplot, 'YLim', [plotmin plotmax]);
    hold(handles.mainplot,'off');
elseif (get(handles.Density, 'Value') == 1)
    % Make density plots
    
    namp = length(snpidx)./5;
    nt = size(snips,2);
        
    % Calculate the mean and sd for each snippet class
    meanSnip = mean(snips(snpidx,:),1);
    sdSnip = std(snips(snpidx,:));
    % Find bounds for histograms
    minAmp = min(meanSnip - 3*sdSnip);
    maxAmp = max(meanSnip + 3*sdSnip);
    % x and y labels for the snip plots
    t=1:nt;
    a=linspace(minAmp, maxAmp, namp);
    
    % make a histogram
    histval = zeros(namp,nt);
    for isn=1:length(snpidx)
        for it=1:nt
            ampind = floor(namp*(snips(snpidx(isn),it)-minAmp)/(maxAmp-minAmp))+1;
            if (ampind <= namp && ampind >= 1)
                histval(ampind,it) = histval(ampind,it) + 1;
            end
        end
    end
    histval = histval./length(snpidx);
    
    pcolor(handles.mainplot, t, a, (histval));
    shading(handles.mainplot, 'interp');
    caxis(handles.mainplot, [0 max(max(histval))/5]);
    colormap('gray');
    axis xy;
    
    hold(handles.mainplot,'on');
    plot(handles.mainplot, meanSnip,'r', 'linewidth', 2);
    plot(handles.mainplot, meanSnip + sdSnip, 'r--');
    plot(handles.mainplot, meanSnip - sdSnip, 'r--');
    hold(handles.mainplot, 'off');
    set(handles.mainplot, 'YLim', [plotmin plotmax]);
else
    for kk=1:kms
        dsnips = bsxfun(@minus, snips(snpidx(labels==kk),:), msnips);
        ktmpv = (r5pc(:,1:2)'*dsnips');
        plot(handles.mainplot, ktmpv(1,:), ktmpv(2,:), '+', 'Color', colorvals(kk));
        hold(handles.mainplot,'on');
    end
    hold(handles.mainplot,'off');
end

%Now draw the individual graphs
    for kk = 1:kms
        hh = ['handles.kmean' num2str(kk)];
        bdf = get(eval(hh), 'ButtonDownFcn');
        cla(eval(hh), 'reset')
        plot(eval(hh), snips(snpidx(labels==kk),:)', 'Color', colorvals(kk));
        set(eval(hh), 'YLim', [plotmin plotmax]);
        set(eval(hh), 'ButtonDownFcn', bdf)
    end

% --- Executes on button press in undo.
function undo_Callback(hObject, eventdata, handles)
% hObject    handle to undo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels validtrials tmpv tmpv2 tmpclusters r5pc s5pc 
global heightsmean heightsstd widthsmean widthsstd goodflag fname tmpsortidx tmpunsortedidx clusize vs tridx tmps2 retry undohistory snips shuffle

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
tempcluster_tlab
set(handles.figure1, 'Visible', 'on')

kms = 6;


snpidx = intersect(shuffle(snpidx1:snpidx2), tmpunsortedidx);

if ~isempty(snpidx)

    if get(handles.slope, 'Value')==1
        labels = cluster(diff(snips(snpidx,:)')', kms, handles);
    else
        labels = cluster(snips(snpidx,:), kms, handles);
    end
    plot_main(handles);
 
end


% --- Executes on button press in selectfile.
function selectfile_Callback(hObject, eventdata, handles)
% hObject    handle to selectfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global unit snips fname nclasses nfiles shuffle


ntsnip = 18;       % Lenth of our spike snippets

[FileName,PathName] = uigetfile('*.h5','Select the hdf5 code file');
bdf2 = disable_all(handles);
fname = [PathName FileName];

fprintf(1,'Reading data from h5 file %s\n', fname);
unit = read_unit_h5file(fname, 'r');
fprintf(1,'...Done\n');

source_directory = unit.source_directory;
% The source directories are in different location on the local computers.
if ismac()
    [status username] = system('who am i');
    if strcmp(strtok(username), 'frederictheunissen')
        if strncmp('/auto/fdata/solveig',source_directory, 19)
            source_directory = strcat('/Users/frederictheunissen/Documents/Data/solveig', source_directory(20:end));
        elseif strncmp('/auto/fdata/julie',source_directory, 17)
            source_directory = strcat('/Users/frederictheunissen/Documents/Data/Julie', source_directory(18:end));
        end
    elseif strcmp(strtok(username), 'elie')
        if strncmp('/auto/fdata/solveig',source_directory, 19)
            source_directory = strcat('/Users/frederictheunissen/Documents/Data/solveig', source_directory(20:end));
        elseif strncmp('/auto/fdata/julie',source_directory, 17)
            source_directory = strcat('/Users/elie/Documents/MATLAB/data', source_directory(18:end));
        end
    elseif strcmp(strtok(username), 'Solveig')
        if strncmp('/auto/fdata/solveig',source_directory, 19)
            source_directory = strcat('/Users/Solveig/PhD', source_directory(20:end));
        elseif strncmp('/auto/k8/fdata/solveig/' ,source_directory, 23)
            source_directory = strcat('/Users/Solveig/PhD/ELECTROPHY/Export_Data', source_directory(34:end));
        elseif strncmp('/auto/fdata/julie',source_directory, 17)
            source_directory = strcat('/Users/frederictheunissen/Documents/Data/Julie', source_directory(18:end));
        end
    end
end


% Loop through data structure and get some stats
nfilesTot = 0;
ntrialsTot = 0;
nsnipsTot =0;

nclasses = length(unit.classes);
nfiles = zeros(nclasses,1);
for iclass=1:nclasses
    nfiles(iclass) = length(unit.class_responses.(unit.classes{iclass})); 
    nfilesTot = nfilesTot + nfiles(iclass);     
    for nfi=1:nfiles(iclass)   
        response = unit.class_responses.(unit.classes{iclass}){nfi};    
        ntrials = length(response.trials); 
        ntrialsTot = ntrialsTot + ntrials;
        
        for it=1:ntrials
            trial = response.trials{it};
            spike_id_trials = trial.spikeIds;
            ns = length(spike_id_trials);
            nsnipsTot = nsnipsTot + ns;
        end
    end
end
fprintf(1, 'File statistics:\n\tNumber of Stims %d\n\tNumber of trials %d\n\tNumber of snips %d\n', nfilesTot, ntrialsTot, nsnipsTot);

% Read all the snippets

snips = zeros(nsnipsTot,ntsnip);

cnt=0;
for iclass=1:nclasses
    a = dir([source_directory '/' unit.site '*' unit.classes{iclass} '*waves.f32']);
    fprintf(1, 'Reading and storing snipets from wave file: %s\n',  [source_directory '/' a.name]);
    fsnip = fopen([source_directory '/' a.name], 'r');
    for nfi=1:nfiles(iclass)   
        response = unit.class_responses.(unit.classes{iclass}){nfi};    
        ntrials = length(response.trials); 
        
        for it=1:ntrials
            trial = response.trials{it};
            spike_id_trials = trial.spikeIds;
            ns = length(spike_id_trials);
            for is=1:ns
                cnt = cnt+1;
                fseek(fsnip, (spike_id_trials(is)-1)*4*ntsnip, 'bof');
                snips(cnt,:) = fread(fsnip, ntsnip, 'float32');
                % Note that we must use the full name if we want to store
                % it.
                unit.class_responses.(unit.classes{iclass}){nfi}.trials{it}.spikeWav(is,:) = snips(cnt,:);
                if mod(cnt, 5000)==0
                    fprintf('.')
                end
            end
        end
    end
    fclose(fsnip);
end

shuffle = randperm(size(snips,1));

fprintf(1, '\n...Done reading and storing snippets\n');

enable_all(handles, bdf2);
sort_function(handles);


% --- Executes on button press in sortall.
function sortall_Callback(hObject, eventdata, handles)
% hObject    handle to sortall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global retry undohistory
bdf2 = disable_all(handles);
fc = finalcluster_tlab;
uiwait(fc)
enable_all(handles, bdf2)
if retry
    clear undohistory
    sort_function(handles)
    tempcluster_tlab
else
    close all
end
    


% --- Executes on button press in Traces.
function Traces_Callback(hObject, eventdata, handles)
% hObject    handle to Traces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Traces
set(hObject, 'Value', 1);
set(handles.Density, 'Value',0);
set(handles.PCA, 'Value',0);
plot_main(handles);

% --- Executes on button press in Density.
function Density_Callback(hObject, eventdata, handles)
% hObject    handle to Density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Density
set(hObject, 'Value', 1);
set(handles.Traces, 'Value',0);
set(handles.PCA, 'Value',0);
plot_main(handles);


% --- Executes on button press in PCA.
function PCA_Callback(hObject, eventdata, handles)
% hObject    handle to PCA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PCA
set(hObject, 'Value', 1);
set(handles.Traces, 'Value',0);
set(handles.Density, 'Value',0);
plot_main(handles);
