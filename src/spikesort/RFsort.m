function varargout = RFsort(varargin)
% RFSORT MATLAB code for RFsort.fig
%      RFSORT, by itself, creates a new RFSORT or raises the existing
%      singleton*.
%
%      H = RFSORT returns the handle to a new RFSORT or the handle to
%      the existing singleton*.
%
%      RFSORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RFSORT.M with the given input arguments.
%
%      RFSORT('Property','Value',...) creates a new RFSORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RFsort_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RFsort_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RFsort

% Last Modified by GUIDE v2.5 25-Apr-2012 10:59:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RFsort_OpeningFcn, ...
                   'gui_OutputFcn',  @RFsort_OutputFcn, ...
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


% --- Executes just before RFsort is made visible.
function RFsort_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RFsort (see VARARGIN)
global sortlisthandle batchnum
% Choose default command line output for RFsort
handles.output = hObject;
batchnum = 0;
% Update handles structure
guidata(hObject, handles);
% plotall(handles)
sortlisthandle = handles.sortlist;

javaaddpath('/usr/share/java/mysql-connector-java.jar');
addpath('/auto/k1/shinji/matlab/arraydb')



% fname = get(handles.filename, 'String');



function plotall(handles)
global batchnum sortlisthandle
% for ii=1:16
%     plot(eval(['handles.sort' num2str(ii)]), randn(10,100));
% end
global partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv goodflag tmpclusters fname sortparams
global heightsmean heightsstd widthsmean widthsstd s5pc r5pc ens tmpsortidx clusize tmpv2 tmpunsortedidx vs tridx tmps2 retry undohistory colortable resorter_handles
global plotmin plotmax
sortlisthandle = handles.sortlist;
% bdf = get(handles.allclustsort, 'ButtonDownFcn');

% cla(handles.allclust, 'reset')
% cla(handles.allclustsort, 'reset')
% hold(handles.allclust);

for cc = 1:16
    ch = cc+batchnum;
    
    if ch > 96
        delete(handles.figure1);
        RFresorter
        return
    end
    
    tmpv2 = zeros(size(cat(2,vv{uvsorts, ch, :})));
    tmps2 = zeros(1,size(tmpv2,2));
    tridx = [];
    uvidx = 0;
    clusize = [];
    uvsize=[];
    
    [vsorts, vtrials] = find(cellfun(@(x) ~isempty(x),vv(:,ch,:))==1);

    uvsorts = unique(vsorts);
    

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

    if length(clusize)>16
        clusize = clusize(1:16);
    end


    

    tmpv3 = [];
    allrange=[];
    
    cla(eval(['handles.sort' num2str(cc)]), 'reset')
    hold(eval(['handles.sort' num2str(cc)]));
    
    colors = sort(randsample(16, length(clusize)));
    colors = colortable(colors);
    
    for ii = 1:length(clusize)
        if ii == 1
            range = 1:clusize(ii);
        else
            range = clusize(ii-1)+1:clusize(ii);
        end

        range = randsample(range, min([length(range), 50]));

        tmpv3 = [tmpv3 tmpv2(:,range)];
        allrange = [allrange range];
        
        plot(eval(['handles.sort' num2str(cc)]), tmpv2(:,range), 'color', colors{ii}./255); set(eval(['handles.sort' num2str(cc)]), 'YLim', [plotmin plotmax]);
    %     hold(handles.allclust); 
        plot(eval(['handles.sort' num2str(cc)]), mean(tmpv2(:,range),2), 'k', 'linewidth', 5)
        plot(eval(['handles.sort' num2str(cc)]), mean(tmpv2(:,range),2), 'color', colors{ii}./255, 'linewidth', 1)
    end
    drawnow
end



% UIWAIT makes RFsort wait for user response (see UIRESUME)
% uiwait(handles.figure1);
for ii=1:16
    set(eval(['handles.sort' num2str(ii)]), 'UserData', 0);
    set(eval(['handles.sort' num2str(ii)]), 'ButtonDownFcn', @select);
    set(eval(['handles.sort' num2str(ii)]), 'Tag', num2str(batchnum + ii));
end

% --- Outputs from this function are returned to the command line.
function varargout = RFsort_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in nextset.
function nextset_Callback(hObject, eventdata, handles)
% hObject    handle to nextset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global batchnum
batchnum = batchnum + 16;
plotall(handles)


function select(src,event)
global selected sortlisthandle

chnum = str2num(get(src, 'Tag'));

if get(src, 'UserData') == 0
    set(src, 'UserData', 1);
    set(src, 'Color', [1 .9 .9])
    selected = sort([selected chnum]);
else
    set(src, 'UserData', 0);
    set(src, 'Color', [1 1 1])
    selected = setdiff(selected, chnum);
end

set(sortlisthandle, 'String', num2str(selected))



function filenum_Callback(hObject, eventdata, handles)
% hObject    handle to filenum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filenum as text
%        str2double(get(hObject,'String')) returns contents of filenum as a double


% --- Executes during object creation, after setting all properties.
function filenum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filenum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadbutton.
function loadbutton_Callback(hObject, eventdata, handles)
% hObject    handle to loadbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sortlisthandle batchnum filenum fname

basedir='/auto/k6/shinji/ziggy_chronic01/';

d=arraydbgettask('runclass','FP30');
batchnum = 0;

% fname = get(handles.filename, 'String');


filenum=str2num(get(handles.filenum, 'String'));
fname=[d(filenum).rawtaskpath d(filenum).rawtaskfile '.tdt.mat'];

set(handles.filename, 'String', fname)

if exist(fname) ~= 2
    error('File does not exist')
    return
end

s1 = strfind(fname, '/');

if exist([fname(1:s1(end)) 'sortparams.mat'])
    error('Already Sorted')
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

global partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv goodflag tmpclusters sortparams
global heightsmean heightsstd widthsmean widthsstd s5pc r5pc ens tmpsortidx clusize tmpv2 tmpunsortedidx vs tridx tmps2 retry undohistory colortable resorter_handles
global plotmin plotmax
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

[plotmin, plotmax] = findplotminmaxstd(cat(2,vv{:, :, :}));

plotall(handles)
