function varargout = tempcluster(varargin)
% TEMPCLUSTER MATLAB code for tempcluster.fig
%      TEMPCLUSTER, by itself, creates a new TEMPCLUSTER or raises the existing
%      singleton*.
%
%      H = TEMPCLUSTER returns the handle to a new TEMPCLUSTER or the handle to
%      the existing singleton*.
%
%      TEMPCLUSTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TEMPCLUSTER.M with the given input arguments.
%
%      TEMPCLUSTER('Property','Value',...) creates a new TEMPCLUSTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tempcluster_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tempcluster_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tempcluster

% Last Modified by GUIDE v2.5 09-Mar-2012 10:26:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tempcluster_OpeningFcn, ...
                   'gui_OutputFcn',  @tempcluster_OutputFcn, ...
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


% --- Executes just before tempcluster is made visible.
function tempcluster_OpeningFcn(hObject, eventdata, handles, varargin)
global partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv2 tmpsortidx
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tempcluster (see VARARGIN)

% Choose default command line output for tempcluster
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% tmpv2(:,snpidx(labels==kk))

% plts = unique(tmpsortidx(tmpsortidx>0));

% for vt = validtrials
%     szs(vt) = size(tmpclusters(vt).v,2);
% end

% ss = find(szs>=ch);
% 
% if ~isempty(ss)

%     tcv = cat(3, tmpclusters(ss).v);

    % tmpcv = cat(2,tcv{:, ch, :});

    for kk = 1:6
        hh = ['handles.cluster' num2str(kk)];
    %     bdf = get(eval(hh), 'ButtonDownFcn');
        cla(eval(hh), 'reset')
        plot(eval(hh), tmpv2(:,tmpsortidx==kk)); set(eval(hh), 'YLim', [plotmin plotmax]);
        hold(eval(hh)); plot(eval(hh), mean(tmpv2(:,tmpsortidx==kk),2), 'k', 'linewidth', 4);
    %     set(eval(hh), 'ButtonDownFcn', bdf)
    end
% else
%     for kk = 1:6
%         hh = ['handles.cluster' num2str(kk)];
%     %     bdf = get(eval(hh), 'ButtonDownFcn');
%         cla(eval(hh), 'reset')
%         plot(eval(hh), []); set(eval(hh), 'YLim', [plotmin plotmax]);
%     %     set(eval(hh), 'ButtonDownFcn', bdf)
%     end
% end





% UIWAIT makes tempcluster wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tempcluster_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
