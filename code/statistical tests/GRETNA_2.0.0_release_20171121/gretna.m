function varargout = gretna(varargin)
% GRETNA MATLAB code for gretna.fig
%      GRETNA, by itself, creates a new GRETNA or raises the existing
%      singleton*.
%
%      H = GRETNA returns the handle to a new GRETNA or the handle to
%      the existing singleton*.
%
%      GRETNA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GRETNA.M with the given input arguments.
%
%      GRETNA('Property','Value',...) creates a new GRETNA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gretna_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gretna_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gretna

% Last Modified by GUIDE v2.5 27-Feb-2017 22:55:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gretna_OpeningFcn, ...
                   'gui_OutputFcn',  @gretna_OutputFcn, ...
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


% --- Executes just before gretna is made visible.
function gretna_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gretna (see VARARGIN)
GRETNAPath=fileparts(which('gretna.m'));
vfid=fopen(fullfile(GRETNAPath, 'VERSION'), 'r');
v=textscan(vfid, '%s');
fclose(vfid);
VERSION=v{1}{1};

fprintf('\tThank you for using GRETNA (v%s)\n\n', VERSION);
fprintf('\tPlease cite:\n') 
fprintf('\tWang, J., Wang, X., Xia, M., Liao, X., Evans, A., & He, Y. (2015).\n')
fprintf('\t\tGRETNA: a graph theoretical network analysis toolbox for imaging connectomics.\n')
fprintf('\t\tFrontiers in human neuroscience, 9, 386.\n');
% Choose default command line output for gretna
handles.output = hObject;
axes(handles.LogoAxes);
axis image;
GRETNAPath=fileparts(which('gretna.m'));
[A, map]=imread(fullfile(GRETNAPath, 'logo.png'));
h=imshow(A, map);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gretna wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gretna_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in PreprocessButton.
function PreprocessButton_Callback(hObject, eventdata, handles)
% hObject    handle to PreprocessButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gretna_GUI_PreprocessInterface;

% --- Executes on button press in CalButton.
function CalButton_Callback(hObject, eventdata, handles)
% hObject    handle to CalButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gretna_GUI_NetAnalysisInterface;


% --- Executes on button press in CompButton.
function CompButton_Callback(hObject, eventdata, handles)
% hObject    handle to CompButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gretna_GUI_CompInterface;


% --- Executes on button press in PlotButton.
function PlotButton_Callback(hObject, eventdata, handles)
% hObject    handle to PlotButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gretna_GUI_PlotBBVInterface;


% --- Executes on button press in GannmBtn.
function GannmBtn_Callback(hObject, eventdata, handles)
% hObject    handle to GannmBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gretna_GUI_GannmInterface;
