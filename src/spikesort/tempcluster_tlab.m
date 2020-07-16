function varargout = tempcluster_tlab(varargin)
% TEMPCLUSTER_TLAB MATLAB code for tempcluster_tlab.fig
%      TEMPCLUSTER_TLAB, by itself, creates a new TEMPCLUSTER_TLAB or raises the existing
%      singleton*.
%
%      H = TEMPCLUSTER_TLAB returns the handle to a new TEMPCLUSTER_TLAB or the handle to
%      the existing singleton*.
%
%      TEMPCLUSTER_TLAB('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TEMPCLUSTER_TLAB.M with the given input arguments.
%
%      TEMPCLUSTER_TLAB('Property','Value',...) creates a new TEMPCLUSTER_TLAB or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tempcluster_tlab_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tempcluster_tlab_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tempcluster_tlab

% Last Modified by GUIDE v2.5 28-Jun-2012 17:25:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tempcluster_tlab_OpeningFcn, ...
                   'gui_OutputFcn',  @tempcluster_tlab_OutputFcn, ...
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


% --- Executes just before tempcluster_tlab is made visible.
function tempcluster_tlab_OpeningFcn(hObject, eventdata, handles, varargin)
global partialflag snpidx snpidx1 snpidx2 plotmin plotmax PF vv ch currsort uvsorts vtrials vsorts labels csort validtrials tmpv2 tmpsortidx snips
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tempcluster_tlab (see VARARGIN)

% Choose default command line output for tempcluster_tlab
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



    for kk = 1:6
        hh = ['handles.cluster' num2str(kk)];
        cla(eval(hh), 'reset')
        plot(eval(hh), snips(tmpsortidx==kk,:)'); set(eval(hh), 'YLim', [plotmin plotmax]);
        hold(eval(hh)); plot(eval(hh), mean(snips(tmpsortidx==kk,:)',2), 'k', 'linewidth', 4);

    end


% UIWAIT makes tempcluster_tlab wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tempcluster_tlab_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
