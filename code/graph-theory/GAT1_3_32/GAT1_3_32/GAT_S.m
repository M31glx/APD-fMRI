function varargout = GAT_S(varargin)
% GAT_S MATLAB code for GAT_S.fig
%      GAT_S, by itself, creates a new GAT_S or raises the existing
%      singleton*.
%
%      H = GAT_S returns the handle to a new GAT_S or the handle to
%      the existing singleton*.
%
%      GAT_S('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GAT_S.M with the given input arguments.
%
%      GAT_S('Property','Value',...) creates a new GAT_S or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GAT_S_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GAT_S_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GAT_S

% Last Modified by GUIDE v2.5 12-Jan-2012 16:56:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GAT_S_OpeningFcn, ...
                   'gui_OutputFcn',  @GAT_S_OutputFcn, ...
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


% --- Executes just before GAT_S is made visible.
function GAT_S_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GAT_S (see VARARGIN)

% Choose default command line output for GAT_S
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GAT_S wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GAT_S_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in ExtractROIs.
function ExtractROIs_Callback(hObject, eventdata, handles)
% hObject    handle to ExtractROIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%global GroupName
%GroupName=input('type name of the group in single quotation: ');
rex;
%rex_out=load('REX.mat');%output from REX
%mat4GAT=params.ROIdata;
%regionName=params.ROInames;
%save(['mat4GAT_' GroupName],'mat4GAT','regionName');

% --- Executes on button press in LoadData.
function LoadData_Callback(hObject, eventdata, handles)
% hObject    handle to LoadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Load_data


% --- Executes on button press in GraphMes_MinDens.
function GraphMes_MinDens_Callback(hObject, eventdata, handles)
% hObject    handle to GraphMes_MinDens (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
KAN_FixedDensity


% --- Executes on button press in GenRandomGraphs.
function GenRandomGraphs_Callback(hObject, eventdata, handles)
% hObject    handle to GenRandomGraphs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Rand_Net_Bootstrap


% --- Executes on button press in CompGraph_AcrossDens.
function CompGraph_AcrossDens_Callback(hObject, eventdata, handles)
% hObject    handle to CompGraph_AcrossDens (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
NetMes_vs_Density_vsNull_final_binary_Corrected_GUI


% --- Executes on button press in CompGraphs_MinDens.
function CompGraphs_MinDens_Callback(hObject, eventdata, handles)
% hObject    handle to CompGraphs_MinDens (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PvalueAtExactMinThreshold


% --- Executes on button press in RegionalMes.
function RegionalMes_Callback(hObject, eventdata, handles)
% hObject    handle to RegionalMes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Net_RegDiff


% --- Executes on button press in NetHubs.
function NetHubs_Callback(hObject, eventdata, handles)
% hObject    handle to NetHubs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Net_Hubs


% --- Executes on button press in TargAttack.
function TargAttack_Callback(hObject, eventdata, handles)
% hObject    handle to TargAttack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
TargAttack_AllNet_AllNetMeasures


% --- Executes on button press in RandFailure.
function RandFailure_Callback(hObject, eventdata, handles)
% hObject    handle to RandFailure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RandAttack_AllNet


% --- Executes on button press in NetModules.
function NetModules_Callback(hObject, eventdata, handles)
% hObject    handle to NetModules (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Module_Analysis


% --- Executes on button press in CompModules.
function CompModules_Callback(hObject, eventdata, handles)
% hObject    handle to CompModules (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Compare_Modules


% --- Executes on button press in NetHubVis.
function NetHubVis_Callback(hObject, eventdata, handles)
% hObject    handle to NetHubVis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
NetHubVisualization


% --- Executes on button press in RegMaps.
function RegMaps_Callback(hObject, eventdata, handles)
% hObject    handle to RegMaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

RegionalMapping


% --- Executes on button press in DegDist.
function DegDist_Callback(hObject, eventdata, handles)
% hObject    handle to DegDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Deg_Distribution

% --- Executes on button press in IndivCont.
function IndivCont_Callback(hObject, eventdata, handles)
% hObject    handle to IndivCont (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Find_IndivContribution 


% --- Executes on button press in AUC_NetMes.
function AUC_NetMes_Callback(hObject, eventdata, handles)
% hObject    handle to AUC_NetMes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AUC_Switch

% --- Executes on button press in AUC_AttackAna.
function AUC_AttackAna_Callback(hObject, eventdata, handles)
% hObject    handle to AUC_AttackAna (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

AUCvsFDA = input('Type 1 for AUC or 2 for FDA analysis: ');

switch AUCvsFDA
    
    case 1
        
        AUC_Analysis_Attack
        
    case 2
        
        FDA_Analysis_Attack
end
    
    
    
