function varargout = GAT_D(varargin)
% GAT_D MATLAB code for GAT_D.fig
%      GAT_D, by itself, creates a new GAT_D or raises the existing
%      singleton*.
%
%      H = GAT_D returns the handle to a new GAT_D or the handle to
%      the existing singleton*.
%
%      GAT_D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GAT_D.M with the given input arguments.
%
%      GAT_D('Property','Value',...) creates a new GAT_D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GAT_D_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GAT_D_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GAT_D

% Last Modified by GUIDE v2.5 28-Jan-2013 16:43:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GAT_D_OpeningFcn, ...
                   'gui_OutputFcn',  @GAT_D_OutputFcn, ...
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


% --- Executes just before GAT_D is made visible.
function GAT_D_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GAT_D (see VARARGIN)

% Choose default command line output for GAT_D
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GAT_D wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GAT_D_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in NetMesf.
function NetMesf_Callback(hObject, eventdata, handles)
% hObject    handle to NetMesf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
NetMes_vs_Density_vsNull_final_binary_Corrected_GUI_D


% --- Executes on button press in CompareNets.
function CompareNets_Callback(hObject, eventdata, handles)
% hObject    handle to CompareNets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Rand_Net_Bootstrap_AcrossSbjComp_d



% --- Executes on button press in TestDisc.
function TestDisc_Callback(hObject, eventdata, handles)
% hObject    handle to TestDisc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Test_Discon_NetMesBin_d

% --- Executes on button press in CompAvg.
function CompAvg_Callback(hObject, eventdata, handles)
% hObject    handle to CompAvg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
average_across_densities_d


% --- Executes on button press in AUC_NetMesf.
function AUC_NetMesf_Callback(hObject, eventdata, handles)
% hObject    handle to AUC_NetMesf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

AUCvsFDA = input('Type 1 for AUC or 2 for FDA analysis: ');

switch AUCvsFDA
    
    case 1
        
        AUC_Analysis_NetMesd
        
    case 2
        
        FDA_Analysis_NetMesd
end



% --- Executes on button press in RegMesf.
function RegMesf_Callback(hObject, eventdata, handles)
% hObject    handle to RegMesf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Net_RegDiff_f_withBootstrap_d
%the above function also performs AUC analysis


% --- Executes on button press in NetHubsf.
function NetHubsf_Callback(hObject, eventdata, handles)
% hObject    handle to NetHubsf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Analysis_Type = input ('hubs based on AUC (type 1) or  FDA (type 2)? ');

switch Analysis_Type
    
    case 1%AUC
        
        AUC_Net_Hubs_d
                
    case 2%FDA
        
        FDA_Net_Hubs_d
        
end




% --- Executes on button press in DegDist.
function DegDist_Callback(hObject, eventdata, handles)
% hObject    handle to DegDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Deg_Distribution_d



% --- Executes on button press in NetHubsVis.
function NetHubsVis_Callback(hObject, eventdata, handles)
% hObject    handle to NetHubsVis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

NetHubVisualization_d


% --- Executes on button press in RegMapf.
function RegMapf_Callback(hObject, eventdata, handles)
% hObject    handle to RegMapf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Regional_Map_d

% --- Executes on button press in NetMod.
function NetMod_Callback(hObject, eventdata, handles)
% hObject    handle to NetMod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Module_Analysis_d 


% --- Executes on button press in CompMod.
function CompMod_Callback(hObject, eventdata, handles)
% hObject    handle to CompMod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Compare_Modules_d