function varargout = GAT(varargin)
% GAT MATLAB code for GAT.fig
%      GAT, by itself, creates a new GAT or raises the existing
%      singleton*.
%
%      H = GAT returns the handle to a new GAT or the handle to
%      the existing singleton*.
%
%      GAT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GAT.M with the given input arguments.
%
%      GAT('Property','Value',...) creates a new GAT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GAT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GAT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GAT

% Last Modified by GUIDE v2.5 15-Dec-2011 16:42:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GAT_OpeningFcn, ...
                   'gui_OutputFcn',  @GAT_OutputFcn, ...
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


% --- Executes just before GAT is made visible.
function GAT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GAT (see VARARGIN)

% Choose default command line output for GAT
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GAT wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GAT_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in VBM.
function VBM_Callback(hObject, eventdata, handles)
% hObject    handle to VBM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(gcf)
GAT_S


% --- Executes on button press in fMRI.
function fMRI_Callback(hObject, eventdata, handles)
% hObject    handle to fMRI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(gcf)
GAT_F

% --- Executes on button press in DWI.
function DWI_Callback(hObject, eventdata, handles)
% hObject    handle to DWI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(gcf)
GAT_D

% --- Executes on button press in Bx.
function Bx_Callback(hObject, eventdata, handles)
% hObject    handle to Bx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(gcf)
GAT_B

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over VBM.


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over VBM.
%function VBM_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to VBM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
