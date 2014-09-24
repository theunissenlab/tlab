function varargout = allclusters(varargin)
% ALLCLUSTERS MATLAB code for allclusters.fig
%      ALLCLUSTERS, by itself, creates a new ALLCLUSTERS or raises the existing
%      singleton*.
%
%      H = ALLCLUSTERS returns the handle to a new ALLCLUSTERS or the handle to
%      the existing singleton*.
%
%      ALLCLUSTERS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ALLCLUSTERS.M with the given input arguments.
%
%      ALLCLUSTERS('Property','Value',...) creates a new ALLCLUSTERS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before allclusters_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to allclusters_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help allclusters

% Last Modified by GUIDE v2.5 24-Apr-2012 10:24:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @allclusters_OpeningFcn, ...
                   'gui_OutputFcn',  @allclusters_OutputFcn, ...
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


% --- Executes just before allclusters is made visible.
function allclusters_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to allclusters (see VARARGIN)

% Choose default command line output for allclusters
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv goodflag tmpclusters fname sortparams
global heightsmean heightsstd widthsmean widthsstd s5pc r5pc ens tmpsortidx clusize tmpv2 tmpunsortedidx vs tridx tmps2 retry undohistory colortable dat clu allrange

colors = sort(randsample(16, min(16, length(clusize))));
colors = colortable(colors);


bdf = get(handles.allclustsort, 'ButtonDownFcn');

cla(handles.allclust, 'reset')
cla(handles.allclustsort, 'reset')
hold(handles.allclust);
tmpv3 = [];
allrange=[];

for ii = 1:min(16, length(clusize))
    if ii == 1
        range = 1:clusize(ii);
    else
        range = clusize(ii-1)+1:clusize(ii);
    end

    range = randsample(range, min([length(range), 50]));
    
    tmpv3 = [tmpv3 tmpv2(:,range)];
    allrange = [allrange range];
    
    plot(handles.allclust, tmpv2(:,range), 'color', colors{ii}./255); set(handles.allclust, 'YLim', [plotmin plotmax]);
%     hold(handles.allclust); 
    plot(handles.allclust, mean(tmpv2(:,range),2), 'k', 'linewidth', 5)
    plot(handles.allclust, mean(tmpv2(:,range),2), 'color', colors{ii}./255, 'linewidth', 1)
end

% set(handles.allclust, 'ButtonDownFcn', bdf)

%% allclustsort
global cellselect
global selectcolors
selectcolors=colormap(jet(8));
selectcolors=selectcolors([1 7 4 6 3 5 2 8],:);
selectcolors = cat(1,[0.6 0.6 0.6], selectcolors);

if isempty(clu)
    clu = 1;
end

% figure(2);clf;
s=size(tmpv,1);
h=plot(handles.allclustsort, tmpv3,'color',selectcolors(1,:));
set(h,'HitTest','off');
xlim([1 s]);
naxis = axis;
hold on;

set(gca, 'ButtonDownFcn', bdf);
dat = [];
dat.snips=tmpv3;
dat.xy=[];
dat.h=h;
dat.sortcode=zeros(1,size(tmpv2,2));
set(gca, 'UserData', dat);




% UIWAIT makes allclusters wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = allclusters_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on mouse press over axes background.
function allclust_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to allclust (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% global clu
% 
% global partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv goodflag tmpclusters fname sortparams
% global heightsmean heightsstd widthsmean widthsstd s5pc r5pc ens tmpsortidx clusize tmpv2 tmpunsortedidx vs tridx tmps2 retry undohistory colortable dat
% 
% 
% selectcolors=colormap(jet(8));
% selectcolors=selectcolors([1 7 2 5 3 6 4 8],:);
% selectcolors = cat(1,[0.6 0.6 0.6], selectcolors);
% 
% pt=get(hObject,'CurrentPoint');
% dat=get(hObject,'UserData');
% 
% x=pt(1,1)
% y=pt(1,2)
% 
% 
% thisxy=[NaN NaN NaN NaN NaN];
% 
% if ~isempty(dat.xy)
%     thisxy=dat.xy(end,:);
%     if ~isnan(thisxy(3))
%         thisxy=[NaN NaN NaN NaN NaN];
%     end
% end
% 
% if isempty(clu)
%     clu=1;
% %     set(handles.cluster1, 'Value',1)
% %     set(handles.cluster2, 'Value',0)
% %     set(handles.cluster3, 'Value',0)
% %     set(handles.cluster4, 'Value',0)
% %     set(handles.cluster5, 'Value',0)
% %     set(handles.cluster6, 'Value',0)
% %     set(handles.deletesc, 'Value',0)
% %     set(handles.examine, 'Value',0)
% 
% end
% 
% if isnan(thisxy(1)) || isnan(clu)
%     thisxy = [x y NaN NaN clu]
%     dat.xy(end+1,:)=thisxy;
%     set(hObject,'UserData',dat);
% else
%     thiscolor = selectcolors(clu+1,:);
%     thisxy(3:4) = [x y];
%     dat.xy(end,:)=thisxy;
%     f=zeros(1,size(tmpv2(:,tmpunsortedidx),2));
% %     s=size(dat.snips,1);
%     s=size(tmpv2,1);
%     'wait'
% %     tic
%     for ii=1:size(tmpv2(:,tmpunsortedidx),2)
%         [xi yi]=polyxpoly(1:s,tmpv2(:,tmpunsortedidx(ii)),thisxy([1 3]),thisxy([2 4]));
%         if ~isempty(xi)
%             f(ii)=1;
%         end
%     end
% %     toc
% %     'wait'
%     
% %     tic
% %     f2 = cellfun(@(ss) ppwrap(1:s,ss,thisxy([1 3]),thisxy([2 4])),num2cell(tmpv2(:,tmpunsortedidx), 1)) ~= 0;
% %     toc
%     'done'
% %     keyboard
% 
%     
% %     dat.sortcode(logical(f))=clu;
%     tmpsortidx(tmpunsortedidx(logical(f))) = clu;
%     plot([thisxy(1) x], [thisxy(2) y],'color',thiscolor,'linewidth',2);
% %     set(dat.h(logical(f)),'color',thiscolor);
%     
%     set(hObject,'UserData',dat);
%     
% end
% 
% 
% % function xout = ppwrap(s,t,xy1,xy2)
% % 
% %     [xout, ~] = polyxpoly(s,t,xy1,xy2);
% %     if isempty(xout)
% %         xout=0;
% %     end
% % xout=xout(1);
% %     


% --- Executes on mouse press over axes background.
function allclustsort_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to allclustsort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv goodflag fname sortparams
global heightsmean heightsstd widthsmean widthsstd s5pc r5pc ens tmpsortidx clusize tmpv2 tmpunsortedidx vs tridx tmps2 retry undohistory colortable dat clu allrange

global cellselect resorter_handles
global selectcolors


pt=get(hObject,'CurrentPoint');

dat=get(hObject,'UserData');

x=pt(1,1);
y=pt(1,2);


thisxy=[NaN NaN NaN NaN NaN];

if ~isempty(dat.xy)
    thisxy=dat.xy(end,:);
    if ~isnan(thisxy(3))
        thisxy=[NaN NaN NaN NaN NaN];
    end
end

% xiv = NaN(1,size(dat.snips,2));
% yiv = NaN(1,size(dat.snips,2));

if isnan(thisxy(1))
    thisxy = [x y NaN NaN clu];
    dat.xy(end+1,:)=thisxy;
    set(hObject,'UserData',dat);
else
    thiscolor = selectcolors(clu+1,:);
    thisxy(3:4) = [x y];
    dat.xy(end,:)=thisxy;
    f=zeros(1,size(dat.snips,2));
    s=size(dat.snips,1);
    
    bdf = disable_all(handles);

    for ii=1:size(dat.snips,2)
        [xi yi]=polyxpoly(1:s,dat.snips(:,ii),thisxy([1 3]),thisxy([2 4]));
        if ~isempty(xi)
            f(ii)=1;
%             xiv(ii) = xi(1);
%             yiv(ii) = yi(1);
        end
    end

%     keyboard
    dat.sortcode(logical(f))=clu;
    tmpsortidx(allrange(logical(f))) = clu;
    tmpunsortedidx = find(tmpsortidx==0);
    plot([thisxy(1) x], [thisxy(2) y],'color',thiscolor,'linewidth',2);
    set(dat.h(logical(f)),'color',thiscolor);
    
    set(hObject,'UserData',dat);
    enable_all(handles, bdf)
    tempcluster
    refresh_resorter(resorter_handles)
    set(hObject, 'Visible', 'on')
end





function bdf = disable_all(handles)
set(handles.status, 'String', 'Wait');  drawnow
bdf = get(handles.allclustsort, 'ButtonDownFcn');
set(handles.allclustsort, 'ButtonDownFcn', []);


function enable_all(handles, bdf)
set(handles.status, 'String', 'Ready');
set(handles.allclustsort, 'ButtonDownFcn', bdf);
