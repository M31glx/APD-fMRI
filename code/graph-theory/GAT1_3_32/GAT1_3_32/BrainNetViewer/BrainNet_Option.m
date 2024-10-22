function varargout = BrainNet_Option(varargin)
%BrainNet Viewer, a graph-based brain network mapping tool, by Mingrui Xia
%Function for control panel
%-----------------------------------------------------------
%	Copyright(c) 2011
%	State Key Laboratory of Cognitive Neuroscience and Learning, Beijing Normal University
%	Written by Mingrui Xia
%	Mail to Author:  <a href="mingruixia@gmail.com">Mingrui Xia</a>
%   Version 1.0;
%   Date 20110531;
%   Last edited 20111027
%-----------------------------------------------------------
%
% BrainNet_Option MATLAB code for BrainNet_Option.fig
%      BrainNet_Option, by itself, creates a new BrainNet_Option or raises the existing
%      singleton*.
%
%      H = BrainNet_Option returns the handle to a new BrainNet_Option or the handle to
%      the existing singleton*.
%
%      BrainNet_Option('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BrainNet_Option.M with the given input arguments.
%
%      BrainNet_Option('Property','Value',...) creates a new BrainNet_Option or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BrainNet_Option_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BrainNet_Option_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BrainNet_Option

% Last Modified by GUIDE v2.5 26-Oct-2011 19:58:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @BrainNet_Option_OpeningFcn, ...
    'gui_OutputFcn',  @BrainNet_Option_OutputFcn, ...
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


% --- Executes just before BrainNet_Option is made visible.
function BrainNet_Option_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BrainNet_Option (see VARARGIN)

% Choose default command line output for BrainNet_Option
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% h_NV=findobj('Tag','NV_fig');
% h_NV=guihandles(h_NV);
% setappdata(handles.EP_fig,'h_NV',h_NV);

movegui(handles.EP_fig,'center');
set(handles.Apply_button,'Enable','off');
Initialization(handles);

function Initialization(handles)
global EC
global surf
global FLAG
switch EC.lot.view
    case 1
        set(handles.LotP_radiobutton,'Value',1);
        set(handles.LotPS_radiobutton,'Enable','on');
        set(handles.LotPA_radiobutton,'Enable','on');
        set(handles.LotPC_radiobutton,'Enable','on');
    case 2
        set(handles.LotF_radiobutton,'Value',1);
        set(handles.LotPS_radiobutton,'Enable','off');
        set(handles.LotPA_radiobutton,'Enable','off');
        set(handles.LotPC_radiobutton,'Enable','off');
    case 3 %%% Added by Mingrui Xia, 20111027, add for medium views(4 views).
        set(handles.LotM_radiobutton, 'Value', 1);
        set(handles.LotPS_radiobutton,'Enable','off');
        set(handles.LotPA_radiobutton,'Enable','off');
        set(handles.LotPC_radiobutton,'Enable','off');
end
switch FLAG.Loadfile
    case {1,3,7,9}
        [t,tl,tr,vl,vr,h1,w1,cut,cuv]=BrainNet('CutMesh', surf);    
        if t<=cut
            set(handles.LotM_radiobutton, 'Enable', 'off');
            if EC.lot.view==3
                EC.lot.view = 2;
                set(handles.LotF_radiobutton,'Value',1);
                set(handles.LotPS_radiobutton,'Enable','off');
                set(handles.LotPA_radiobutton,'Enable','off');
                set(handles.LotPC_radiobutton,'Enable','off');
            end
        end
    case {2,6}
        set(handles.LotM_radiobutton, 'Enable', 'off');
            if EC.lot.view==3
                EC.lot.view = 2;
                set(handles.LotF_radiobutton,'Value',1);
                set(handles.LotPS_radiobutton,'Enable','off');
                set(handles.LotPA_radiobutton,'Enable','off');
                set(handles.LotPC_radiobutton,'Enable','off');
            end    
end
%%% Added END, 20111027
switch EC.lot.view_direction
    case 1
        set(handles.LotPS_radiobutton,'Value',1);
    case 2
        set(handles.LotPA_radiobutton,'Value',1);
    case 3
        set(handles.LotPC_radiobutton,'Value',1);
end
set(handles.BakC_text,'BackgroundColor',EC.bak.color);
if FLAG.Loadfile==1 || FLAG.Loadfile==3 || FLAG.Loadfile==7
    set(handles.MshC_text,'BackgroundColor',EC.msh.color);
    set(handles.MshA_slider, 'Value',EC.msh.alpha);
    set(handles.MshA_edit, 'String',num2str(EC.msh.alpha,'%6.2f'));
else
    set(handles.MshC_text,'Enable','off');
    set(handles.MshA_slider,'Enable','off');
    set(handles.MshA_edit,'Enable','off');
end
if FLAG.Loadfile==2 || FLAG.Loadfile==3 || FLAG.Loadfile==7 || FLAG.Loadfile==6
    switch EC.lbl
        case 1
            set(handles.LblA_radiobutton,'Value',1);
            set(handles.LblT_slider,'Enable','off');
            set(handles.LblT_edit,'Enable','off');
            set(handles.LblT_popupmenu,'Enable','off');
        case 2
            set(handles.LblN_radiobutton,'Value',1);
            set(handles.LblT_slider,'Enable','off');
            set(handles.LblT_edit,'Enable','off');
            set(handles.LblT_popupmenu,'Enable','off');
        case 3
            set(handles.LblT_radiobutton,'Value',1);
            set(handles.LblT_slider,'Enable','on');
            set(handles.LblT_edit,'Enable','on');
            set(handles.LblT_popupmenu,'Enable','on');
    end
    switch EC.lbl_threshold_type
        case 1
            if EC.lbl_threshold>max(surf.sphere(:,5)) || EC.lbl_threshold<min(surf.sphere(:,5))
                EC.lbl_threshold=min(surf.sphere(:,5));
            end
            set(handles.LblT_popupmenu,'Value',1);
            set(handles.LblT_slider,'Max',max(surf.sphere(:,5)));
            set(handles.LblT_slider,'Min',min(surf.sphere(:,5))-0.001);
            set(handles.LblT_slider,'Value',EC.lbl_threshold);
            set(handles.LblT_edit,'String',num2str(EC.lbl_threshold,'%6.2f'));
        case 2
            if EC.lbl_threshold>max(surf.sphere(:,4)) || EC.lbl_threshold<min(surf.sphere(:,4))
                EC.lbl_threshold=min(surf.sphere(:,4));
            end
            set(handles.LblT_popupmenu,'Value',2);
            set(handles.LblT_slider,'Max',max(surf.sphere(:,4)));
            set(handles.LblT_slider,'Min',min(surf.sphere(:,4))-0.001);
            set(handles.LblT_slider,'Value',EC.lbl_threshold);
            set(handles.LblT_edit,'String',num2str(EC.lbl_threshold,'%6.2f'));
    end
    switch EC.nod.draw
        case 1
            set(handles.NodDA_radiobutton,'Value',1);
            set(handles.NodDT_slider,'Enable','off');
            set(handles.NodDT_edit,'Enable','off');
            set(handles.NodDT_popupmenu,'Enable','off');
        case 2
            set(handles.NodDT_radiobutton,'Value',1);
            set(handles.NodDT_slider,'Enable','on');
            set(handles.NodDT_edit,'Enable','on');
            set(handles.NodDT_popupmenu,'Enable','on');
    end
    switch EC.nod.draw_threshold_type
        case 1
            if EC.nod.draw_threshold>max(surf.sphere(:,5)) || EC.nod.draw_threshold<min(surf.sphere(:,5))
                EC.nod.draw_threshold=min(surf.sphere(:,5));
            end
            set(handles.NodDT_popupmenu,'Value',1);
            set(handles.NodDT_slider,'Max',max(surf.sphere(:,5)));
            set(handles.NodDT_slider,'Min',min(surf.sphere(:,5))-0.001);
            set(handles.NodDT_slider,'Value',EC.nod.draw_threshold);
            set(handles.NodDT_edit,'String',num2str(EC.nod.draw_threshold,'%6.2f'));
        case 2
            if EC.nod.draw_threshold>max(surf.sphere(:,4)) || EC.nod.draw_threshold<min(surf.sphere(:,4))
                EC.nod.draw_threshold=min(surf.sphere(:,4));
            end
            set(handles.NodDT_popupmenu,'Value',2);
            set(handles.NodDT_slider,'Max',max(surf.sphere(:,4)));
            set(handles.NodDT_slider,'Min',min(surf.sphere(:,4))-0.001);
            set(handles.NodDT_slider,'Value',EC.nod.draw_threshold);
            set(handles.NodDT_edit,'String',num2str(EC.nod.draw_threshold,'%6.2f'));
    end
    switch EC.nod.size
        case 1
            set(handles.NodSS_radiobutton,'Value',1);
            set(handles.NodSS_edit,'Enable','on');
            set(handles.NodSV_popupmenu,'Enable','off');
            set(handles.NodST_slider,'Enable','off');
            set(handles.NodST_edit,'Enable','off');
            set(handles.NodST_checkbox,'Enable','off');
        case 2
            set(handles.NodSV_radiobutton,'Value',1);
            set(handles.NodSS_edit,'Enable','off');
            set(handles.NodSV_popupmenu,'Enable','on');
            set(handles.NodST_slider,'Enable','off');
            set(handles.NodST_edit,'Enable','off');
            set(handles.NodST_checkbox,'Enable','on');
            set(handles.NodST_checkbox,'Value',0);
        case 3
            set(handles.NodSV_radiobutton,'Value',1);
            set(handles.NodST_checkbox,'Value',1);
            set(handles.NodSS_edit,'Enable','off');
            set(handles.NodSV_popupmenu,'Enable','on');
            set(handles.NodST_slider,'Enable','on');
            set(handles.NodST_edit,'Enable','on');
            set(handles.NodST_checkbox,'Enable','on');
    end
    set(handles.NodSS_edit,'String',num2str(EC.nod.size_size,'%6.2f'));
    set(handles.NodSV_popupmenu,'Value',EC.nod.size_value);
    if EC.nod.size_threshold>max(surf.sphere(:,5)) || EC.nod.size_threshold<min(surf.sphere(:,5))
        EC.nod.size_threshold=min(surf.sphere(:,5));
    end
    set(handles.NodST_slider,'Max',max(surf.sphere(:,5)));
    set(handles.NodST_slider,'Min',min(surf.sphere(:,5))-0.001);
    set(handles.NodST_slider,'Value',EC.nod.size_threshold);
    set(handles.NodST_edit,'String',num2str(EC.nod.size_threshold,'%6.2f'));
    set(handles.NodSVR_slider,'Value',EC.nod.size_ratio);
    set(handles.NodSVR_edit,'String',num2str(EC.nod.size_ratio,'%6.2f'));
    switch EC.nod.color
        case 1
            set(handles.NodCS_radiobutton,'Value',1);
            set(handles.NodCC_popupmenu,'Enable','off');
            set(handles.NodCM_pushbutton,'Enable','off');
            set(handles.NodCT_slider,'Enable','off');
            set(handles.NodCT_edit,'Enable','off');
            set(handles.NodCS_text,'BackgroundColor',EC.nod.CM(1,:));
        case 2
            set(handles.NodCC_radiobutton,'Value',1);
            set(handles.NodCC_popupmenu,'Enable','on');
            set(handles.NodCM_pushbutton,'Enable','off');
            set(handles.NodCT_slider,'Enable','off');
            set(handles.NodCT_edit,'Enable','off');
        case 3
            set(handles.NodCM_radiobutton,'Value',1);
            set(handles.NodCC_popupmenu,'Enable','off');
            set(handles.NodCM_pushbutton,'Enable','on');
            set(handles.NodCT_slider,'Enable','off');
            set(handles.NodCT_edit,'Enable','off');
        case 4
            set(handles.NodCT_radiobutton,'Value',1);
            set(handles.NodCC_popupmenu,'Enable','off');
            set(handles.NodCM_pushbutton,'Enable','off');
            set(handles.NodCT_slider,'Enable','on');
            set(handles.NodCT_edit,'Enable','on');
            set(handles.NodCTH_text,'BackgroundColor',EC.nod.CM(1,:));
            set(handles.NodCTL_text,'BackgroundColor',EC.nod.CM(64,:));
    end
    set(handles.NodCC_popupmenu,'Value',EC.nod.color_map);
    if EC.nod.color_threshold>max(surf.sphere(:,4)) || EC.nod.color_threshold<min(surf.sphere(:,4))
        EC.nod.color_threshold=min(surf.sphere(:,4));
    end
    set(handles.NodCT_slider,'Max',max(surf.sphere(:,4)));
    set(handles.NodCT_slider,'Min',min(surf.sphere(:,4))-0.001);
    set(handles.NodCT_slider,'Value',EC.nod.color_threshold);
    set(handles.NodCT_edit,'String',num2str(EC.nod.color_threshold,'%6.2f'));
else
    set(handles.NodDA_radiobutton,'Enable','off');
    set(handles.NodDT_radiobutton,'Enable','off');
    set(handles.NodDT_slider,'Enable','off');
    set(handles.NodDT_edit,'Enable','off');
    set(handles.NodDT_popupmenu,'Enable','off');
    set(handles.NodSS_radiobutton,'Enable','off');
    set(handles.NodSS_edit,'Enable','off');
    set(handles.NodSV_radiobutton,'Enable','off');
    set(handles.NodSV_popupmenu,'Enable','off');
    set(handles.NodST_checkbox,'Enable','off');
    set(handles.NodST_slider,'Enable','off');
    set(handles.NodST_edit,'Enable','off');
    set(handles.NodSVR_slider,'Enable','off');
    set(handles.NodSVR_edit,'Enable','off');
    set(handles.NodCS_radiobutton,'Enable','off');
    set(handles.NodCC_radiobutton,'Enable','off');
    set(handles.NodCC_popupmenu,'Enable','off');
    set(handles.NodCM_radiobutton,'Enable','off');
    set(handles.NodCM_pushbutton,'Enable','off');
    set(handles.NodCT_radiobutton,'Enable','off');
    set(handles.NodCT_slider,'Enable','off');
    set(handles.NodCT_edit,'Enable','off');
    set(handles.LblA_radiobutton,'Enable','off');
    set(handles.LblN_radiobutton,'Enable','off');
    set(handles.LblT_radiobutton,'Enable','off');
    set(handles.LblT_slider,'Enable','off');
    set(handles.LblT_edit,'Enable','off');
    set(handles.LblT_popupmenu,'Enable','off');
    set(handles.LblF_button,'Enable','off');
end
if FLAG.Loadfile==7 || FLAG.Loadfile==6
    switch EC.edg.draw
        case 1
            set(handles.EdgDA_radiobutton,'Value',1);
            set(handles.EdgDT_slider,'Enable','off');
            set(handles.EdgDT_edit,'Enable','off');
            set(handles.EdgDS_slider,'Enable','off');
            set(handles.EdgDS_edit,'Enable','off');
            set(handles.EdgDT_checkbox,'Enable','off');
        case 2
            set(handles.EdgDT_radiobutton,'Value',1);
            set(handles.EdgDT_slider,'Enable','on');
            set(handles.EdgDT_edit,'Enable','on');
            set(handles.EdgDS_slider,'Enable','on');
            set(handles.EdgDS_edit,'Enable','on');
            set(handles.EdgDT_checkbox,'Enable','on');
            set(handles.EdgDT_checkbox,'Value',EC.edg.draw_abs);
    end
    switch get(handles.EdgDT_checkbox,'Value')
        case 0
            if EC.edg.draw_threshold>max(max(surf.net)) || EC.edg.draw_threshold<min(min(surf.net))
                EC.edg.draw_threshold=min(min(surf.net));
            end
            set(handles.EdgDT_slider,'Max',max(max(surf.net)));
            set(handles.EdgDT_slider,'Min',min(min(surf.net))-0.001);
            set(handles.EdgDT_slider,'Value',EC.edg.draw_threshold);
            set(handles.EdgDT_edit,'String',num2str(EC.edg.draw_threshold,'%6.2f'));
            set(handles.EdgDS_slider,'Value',length(find(surf.net>EC.edg.draw_threshold))/(size(surf.net,1)*size(surf.net,2)));
            set(handles.EdgDS_edit,'String',num2str(length(find(surf.net>EC.edg.draw_threshold))/(size(surf.net,1)*size(surf.net,2)),'%6.2f'));
        case 1
            if EC.edg.draw_threshold>max(max(abs(surf.net))) || EC.edg.draw_threshold<min(min(abs(surf.net)))
                EC.edg.draw_threshold=min(min(abs(surf.net)));
            end
            set(handles.EdgDT_slider,'Max',max(max(abs(surf.net))));
            set(handles.EdgDT_slider,'Min',min(min(abs(surf.net)))-0.001);
            set(handles.EdgDT_slider,'Value',EC.edg.draw_threshold);
            set(handles.EdgDT_edit,'String',num2str(EC.edg.draw_threshold,'%6.2f'));
            set(handles.EdgDS_slider,'Value',length(find(abs(surf.net)>EC.edg.draw_threshold))/(size(surf.net,1)*size(surf.net,2)));
            set(handles.EdgDS_edit,'String',num2str(length(find(abs(surf.net)>EC.edg.draw_threshold))/(size(surf.net,1)*size(surf.net,2)),'%6.2f'));
    end
    switch EC.edg.size
        case 1
            set(handles.EdgSS_radiobutton,'Value',1);
            set(handles.EdgSS_edit,'Enable','on');
            set(handles.EdgSV_popupmenu,'Enable','off');
            set(handles.EdgST_slider,'Enable','off');
            set(handles.EdgST_edit,'Enable','off');
            set(handles.EdgST_checkbox,'Enable','off');
        case 2
            set(handles.EdgSV_radiobutton,'Value',1);
            set(handles.EdgSS_edit,'Enable','off');
            set(handles.EdgSV_popupmenu,'Enable','on');
            set(handles.EdgST_checkbox,'Enable','on');
            set(handles.EdgST_slider,'Enable','off');
            set(handles.EdgST_edit,'Enable','off');
        case 3
            set(handles.EdgST_checkbox,'Value',1);
            set(handles.EdgSS_edit,'Enable','off');
            set(handles.EdgSV_popupmenu,'Enable','on');
            set(handles.EdgST_slider,'Enable','on');
            set(handles.EdgST_edit,'Enable','on');
            set(handles.EdgSV_radiobutton,'Value',1);
    end
    set(handles.EdgSS_edit,'String',num2str(EC.edg.size_size,'%6.2f'));
    set(handles.EdgSV_popupmenu,'Value',EC.edg.size_value);
    set(handles.EdgS_checkbox,'Value',EC.edg.size_abs);
    switch EC.edg.size_abs
        case 0
            if EC.edg.size_threshold>max(max(surf.net)) || EC.edg.size_threshold<min(min(surf.net))
                EC.edg.size_threshold=min(min(surf.net));
            end
            set(handles.EdgST_slider,'Max',max(max(surf.net)));
            set(handles.EdgST_slider,'Min',min(min(surf.net))-0.001);
            set(handles.EdgST_slider,'Value',EC.edg.size_threshold);
            set(handles.EdgST_edit,'String',num2str(EC.edg.size_threshold,'%6.2f'));
        case 1
            if EC.edg.size_threshold>max(max(abs(surf.net))) || EC.edg.size_threshold<min(min(abs(surf.net)))
                EC.edg.size_threshold=min(min(abs(surf.net)));
            end
            set(handles.EdgST_slider,'Max',max(max(abs(surf.net))));
            set(handles.EdgST_slider,'Min',min(min(abs(surf.net)))-0.001);
            set(handles.EdgST_slider,'Value',EC.edg.size_threshold);
            set(handles.EdgST_edit,'String',num2str(EC.edg.size_threshold,'%6.2f'));
    end
    set(handles.EdgSRR_slider,'Value',EC.edg.size_ratio);
    set(handles.EdgSRR_edit,'String',num2str(EC.edg.size_ratio,'%6.2f'));
    switch EC.edg.color
        case 1
            set(handles.EdgCS_radiobutton,'Value',1);
            set(handles.EdgCS_text,'BackgroundColor',EC.edg.CM(1,:));
            set(handles.EdgCC_popupmenu,'Enable','off');
            set(handles.EdgCT_slider,'Enable','off');
            set(handles.EdgCT_edit,'Enable','off');
            set(handles.EdgCD_slider,'Enable','off');
            set(handles.EdgCD_edit,'Enable','off');
        case 2
            set(handles.EdgCC_radiobutton,'Value',1);
            set(handles.EdgCC_popupmenu,'Enable','on');
            set(handles.EdgCT_slider,'Enable','off');
            set(handles.EdgCT_edit,'Enable','off');
            set(handles.EdgCD_slider,'Enable','off');
            set(handles.EdgCD_edit,'Enable','off');
        case 3
            set(handles.EdgCT_radiobutton,'Value',1);
            set(handles.EdgCC_popupmenu,'Enable','off');
            set(handles.EdgCT_slider,'Enable','on');
            set(handles.EdgCT_edit,'Enable','on');
            set(handles.EdgCT_slider,'Value',EC.edg.color_threshold);
            set(handles.EdgCT_edit,'String',num2str(EC.edg.color_threshold,'%6.2f'));
            set(handles.EdgCTH_text,'BackgroundColor',EC.edg.CM(1,:));
            set(handles.EdgCTL_text,'BackgroundColor',EC.edg.CM(64,:));
            set(handles.EdgCD_slider,'Enable','off');
            set(handles.EdgCD_edit,'Enable','off');
        case 4
            set(handles.EdgCD_radiobutton,'Value',1);
            set(handles.EdgCC_popupmenu,'Enable','off');
            set(handles.EdgCT_slider,'Enable','off');
            set(handles.EdgCT_edit,'Enable','off');
            set(handles.EdgCD_slider,'Enable','on');
            set(handles.EdgCD_edit,'Enable','on');
            set(handles.EdgCDH_text,'BackgroundColor',EC.edg.CM(1,:));
            set(handles.EdgCDL_text,'BackgroundColor',EC.edg.CM(64,:));
            set(handles.EdgCD_slider,'Value',EC.edg.color_distance);
            set(handles.EdgCD_edit,'String',num2str(EC.edg.color_distance,'%6.2f'));
    end
    set(handles.EdgCC_popupmenu,'Value',EC.edg.color_map);
    set(handles.EdgC_checkbox,'Value',EC.edg.color_abs);
    switch EC.edg.color_abs
        case 0
            if EC.edg.color_threshold>max(max(surf.net)) || EC.edg.color_threshold<min(min(surf.net))
                EC.edg.color_threshold=min(min(surf.net));
            end
            set(handles.EdgCT_slider,'Max',max(max(surf.net)));
            set(handles.EdgCT_slider,'Min',min(min(surf.net))-0.001);
            set(handles.EdgCT_slider,'Value',EC.edg.color_threshold);
            set(handles.EdgCT_edit,'String',num2str(EC.edg.color_threshold,'%6.2f'));
        case 1
            if EC.edg.color_threshold>max(max(abs(surf.net))) || EC.edg.color_threshold<min(min(abs(surf.net)))
                EC.edg.color_threshold=min(min(abs(surf.net)));
            end
            set(handles.EdgCT_slider,'Max',max(max(abs(surf.net))));
            set(handles.EdgCT_slider,'Min',min(min(abs(surf.net)))-0.001);
            set(handles.EdgCT_slider,'Value',EC.edg.color_threshold);
            set(handles.EdgCT_edit,'String',num2str(EC.edg.color_threshold,'%6.2f'));
    end
    
else
    set(handles.EdgDA_radiobutton,'Enable','off');
    set(handles.EdgDT_radiobutton,'Enable','off');
    set(handles.EdgDT_slider,'Enable','off');
    set(handles.EdgDT_edit,'Enable','off');
    set(handles.EdgSS_radiobutton,'Enable','off');
    set(handles.EdgSS_edit,'Enable','off');
    set(handles.EdgSV_radiobutton,'Enable','off');
    set(handles.EdgSV_popupmenu,'Enable','off');
    set(handles.EdgST_checkbox,'Enable','off');
    set(handles.EdgST_slider,'Enable','off');
    set(handles.EdgST_edit,'Enable','off');
    set(handles.EdgSRR_slider,'Enable','off');
    set(handles.EdgSRR_edit,'Enable','off');
    set(handles.EdgCS_radiobutton,'Enable','off');
    set(handles.EdgCC_radiobutton,'Enable','off');
    set(handles.EdgCC_popupmenu,'Enable','off');
    set(handles.EdgCT_radiobutton,'Enable','off');
    set(handles.EdgCT_slider,'Enable','off');
    set(handles.EdgCT_edit,'Enable','off');
    set(handles.EdgCD_slider,'Enable','off');
    set(handles.EdgCD_edit,'Enable','off');
    set(handles.EdgDS_slider,'Enable','off');
    set(handles.EdgDS_edit,'Enable','off');
    set(handles.EdgDT_checkbox,'Enable','off');
    set(handles.EdgS_checkbox,'Enable','off');
end
if FLAG.Loadfile==9
    if FLAG.MAP==2
    set(handles.VolDR_text,'String',['Volume Data Range:  ',num2str(min(min(min(surf.mask))),'%6.2f'),'   ',num2str(max(max(max(surf.mask))),'%6.2f')]);
    else
         set(handles.VolDR_text,'String',['Volume Data Range:  ',num2str(min(surf.T),'%6.2f'),'   ',num2str(max(surf.T),'%6.2f')]);
    end
    switch EC.vol.display
        case 1
            set(handles.VolD_popupmenu,'Value',1);
        case 2
            set(handles.VolD_popupmenu,'Value',2);
            set(handles.VolNR_text,'Enable','off');
            set(handles.VolNRn_edit,'Enable','off');
            set(handles.VolNRx_edit,'Enable','off');
        case 3
            set(handles.VolD_popupmenu,'Value',3);
            set(handles.VolPR_text,'Enable','off');
            set(handles.VolPRn_edit,'Enable','off');
            set(handles.VolPRx_edit,'Enable','off');
    end
    set(handles.VolPRn_edit,'String',num2str(EC.vol.pn,'%6.2f'));
    set(handles.VolPRx_edit,'String',num2str(EC.vol.px,'%6.2f'));
    set(handles.VolNRn_edit,'String',num2str(EC.vol.nn,'%6.2f'));
    set(handles.VolNRx_edit,'String',num2str(EC.vol.nx,'%6.2f'));
    set(handles.VolNCS_text,'BackgroundColor',EC.vol.null);
    set(handles.VolC_popupmenu,'Value',EC.vol.color_map);
else
    set(handles.VolDR_text,'Enable','off');
    set(handles.VolD_text,'Enable','off');
    set(handles.VolD_popupmenu,'Enable','off');
    set(handles.VolNR_text,'Enable','off');
    set(handles.VolNRn_edit,'Enable','off');
    set(handles.VolNRx_edit,'Enable','off');
    set(handles.VolPR_text,'Enable','off');
    set(handles.VolPRn_edit,'Enable','off');
    set(handles.VolPRx_edit,'Enable','off');
    set(handles.VolNC_text,'Enable','off');
    set(handles.VolC_text,'Enable','off');
    set(handles.VolC_popupmenu,'Enable','off');
end
set(handles.ImgPW_edit,'String',num2str(EC.img.width));
set(handles.ImgPH_edit,'String',num2str(EC.img.height));
set(handles.ImgD_edit,'String',num2str(EC.img.dpi));
set(handles.ImgDW_edit,'String',num2str(EC.img.width/EC.img.dpi*2.54,'%6.2f'));
set(handles.ImgDH_edit,'String',num2str(EC.img.height/EC.img.dpi*2.54,'%6.2f'));
set(handles.ImgDW_popupmenu,'Value',1);
set(handles.ImgDH_popupmenu,'Value',1);
set(handles.ImgC_checkbox,'Value',1);
if FLAG.LF==1
    set(handles.Apply_button,'Enable','on');
end

% UIWAIT makes BrainNet_Option wait for user response (see UIRESUME)
% uiwait(handles.EP_fig);


% --- Outputs from this function are returned to the command line.
function varargout = BrainNet_Option_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in EC_list.
function EC_list_Callback(hObject, eventdata, handles)
% hObject    handle to EC_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns EC_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from EC_list
str = get(hObject, 'string');
n = get(hObject, 'value');
switch str{n}
    case 'Background'
        set(handles.Msh_panel, 'visible','off');
        set(handles.Bak_panel, 'visible','on');
        set(handles.Nod_panel, 'visible','off');
        set(handles.Edg_panel, 'visible','off');
        set(handles.Img_panel, 'visible','off');
        set(handles.Lot_panel,'visible','off');
        set(handles.Vol_panel,'Visible','off');
    case 'Surface'
        set(handles.Bak_panel, 'visible','off');
        set(handles.Msh_panel, 'visible','on');
        set(handles.Nod_panel, 'visible','off');
        set(handles.Edg_panel, 'visible','off');
        set(handles.Img_panel, 'visible','off');
        set(handles.Lot_panel,'visible','off');
        set(handles.Vol_panel,'Visible','off');
    case 'Node'
        set(handles.Nod_panel,'visible','on');
        set(handles.Msh_panel, 'visible','off');
        set(handles.Bak_panel, 'visible','off');
        set(handles.Edg_panel, 'visible','off');
        set(handles.Img_panel, 'visible','off');
        set(handles.Lot_panel,'visible','off');
        set(handles.Vol_panel,'Visible','off');
    case 'Edge'
        set(handles.Nod_panel,'visible','off');
        set(handles.Msh_panel, 'visible','off');
        set(handles.Bak_panel, 'visible','off');
        set(handles.Edg_panel, 'visible','on');
        set(handles.Img_panel, 'visible','off');
        set(handles.Lot_panel,'visible','off');
        set(handles.Vol_panel,'Visible','off');
    case 'Image'
        set(handles.Nod_panel,'visible','off');
        set(handles.Msh_panel, 'visible','off');
        set(handles.Bak_panel, 'visible','off');
        set(handles.Edg_panel, 'visible','off');
        set(handles.Img_panel, 'visible','on');
        set(handles.Lot_panel,'visible','off');
        set(handles.Vol_panel,'Visible','off');
    case 'Layout'
        set(handles.Nod_panel,'visible','off');
        set(handles.Msh_panel, 'visible','off');
        set(handles.Bak_panel, 'visible','off');
        set(handles.Edg_panel, 'visible','off');
        set(handles.Img_panel, 'visible','off');
        set(handles.Lot_panel,'visible','on');
        set(handles.Vol_panel,'Visible','off');
    case 'Volume'
        set(handles.Nod_panel,'visible','off');
        set(handles.Msh_panel, 'visible','off');
        set(handles.Bak_panel, 'visible','off');
        set(handles.Edg_panel, 'visible','off');
        set(handles.Img_panel, 'visible','off');
        set(handles.Lot_panel,'visible','off');
        set(handles.Vol_panel,'Visible','on');
end

% --- Executes during object creation, after setting all properties.
function EC_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EC_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Reset_button.
function Reset_button_Callback(hObject, eventdata, handles)
% hObject    handle to Reset_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Initialization(handles);


% --- Executes on button press in OK_button.
function OK_button_Callback(hObject, eventdata, handles)
% hObject    handle to OK_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global FLAG
GetValue(handles);
if FLAG.LF==1 || FLAG.EC_change==1
    FLAG.EC=1;
end
close(findobj('Tag','EP_fig'));


function GetValue(handles)
global EC
global FLAG
FLAG.IsCalledByREST = 0;
if get(handles.LotP_radiobutton,'Value')==1
    EC.lot.view = 1;
elseif get(handles.LotF_radiobutton, 'Value')==1
    EC.lot.view = 2;
else %%% Added by Mingrui Xia, 20111026, add for medium views (4 views).
    EC.lot.view = 3;
end
if get(handles.LotPS_radiobutton,'Value')==1
    EC.lot.view_direction=1;
    FLAG.sagittal=1;
elseif get(handles.LotPA_radiobutton,'Value')==1
    EC.lot.view_direction=2;
    FLAG.axis=1;
else
    EC.lot.view_direction=3;
    FLAG.coronal=1;
end
EC.bak.color=get(handles.BakC_text,'BackgroundColor');
EC.msh.color=get(handles.MshC_text,'BackgroundColor');
EC.msh.alpha=str2double(get(handles.MshA_edit,'String'));
if get(handles.NodDA_radiobutton,'Value')==1
    EC.nod.draw=1;
else
    EC.nod.draw=2;
    EC.nod.draw_threshold_type=get(handles.NodDT_popupmenu,'Value');
    EC.nod.draw_threshold=str2double(get(handles.NodDT_edit,'String'));
end

if get(handles.NodSS_radiobutton,'Value')==1
    EC.nod.size=1;
    EC.nod.size_size=str2double(get(handles.NodSS_edit,'String'));
elseif get(handles.NodSV_radiobutton,'Value')==1
    if get(handles.NodST_checkbox,'Value')==1
        EC.nod.size=3;
        EC.nod.size_threshold=str2double(get(handles.NodST_edit,'String'));
    else
        EC.nod.size=2;
        EC.nod.size_value=get(handles.NodSV_popupmenu,'Value');
    end
end

EC.nod.size_ratio=str2double(get(handles.NodSVR_edit,'String'));
if get(handles.NodCS_radiobutton,'Value')==1
    EC.nod.color=1;
    EC.nod.CM(1,:)=get(handles.NodCS_text,'BackgroundColor');
elseif get(handles.NodCC_radiobutton,'Value')==1
    EC.nod.color=2;
    EC.nod.color_map=get(handles.NodCC_popupmenu,'Value');
    switch EC.nod.color_map
        case 1
            EC.nod.CM=jet;
        case 2
            EC.nod.CM=hsv;
        case 3
            EC.nod.CM=hot;
        case 4
            EC.nod.CM=cool;
        case 5
            EC.nod.CM=spring;
        case 6
            EC.nod.CM=summer;
        case 7
            EC.nod.CM=autumn;
        case 8
            EC.nod.CM=winter;
        case 9
            EC.nod.CM=gray;
        case 10
            EC.nod.CM=bone;
        case 11
            EC.nod.CM=copper;
        case 12
            EC.nod.CM=pink;
        case 13
            EC.nod.CM=lines;
    end
elseif get(handles.NodCM_radiobutton,'Value')==1
    EC.nod.color=3;
    EC.nod.CM(1:21,:)=EC.nod.CMm;
    EC.nod.CM(22,:)=[0.5,0.5,0.5];
else
    EC.nod.color=4;
    EC.nod.color_threshold=str2double(get(handles.NodCT_edit,'String'));
    EC.nod.CM(1,:)=get(handles.NodCTH_text,'BackgroundColor');
    EC.nod.CM(64,:)=get(handles.NodCTL_text,'BackgroundColor');
end
if get(handles.LblA_radiobutton,'Value')==1
    EC.lbl=1;
elseif get(handles.LblN_radiobutton,'Value')==1
    EC.lbl=2;
else
    EC.lbl=3;
    EC.lbl_threshold=str2double(get(handles.LblT_edit,'String'));
    EC.lbl_threshold_type=get(handles.LblT_popupmenu,'Value');
end
if get(handles.EdgDA_radiobutton,'Value')==1
    EC.edg.draw=1;
else
    EC.edg.draw=2;
    EC.edg.draw_threshold=str2double(get(handles.EdgDT_edit,'String'));
    EC.edg.draw_abs=get(handles.EdgDT_checkbox,'Value');
end
if get(handles.EdgSS_radiobutton,'Value')==1
    EC.edg.size=1;
    EC.edg.size_size=str2double(get(handles.EdgSS_edit,'String'));
elseif get(handles.EdgSV_radiobutton,'Value')==1
    if get(handles.EdgST_checkbox,'Value')==1
        EC.edg.size=3;
        EC.edg.size_threshold=str2double(get(handles.EdgST_edit,'String'));
    else
        EC.edg.size=2;
        EC.edg.size_value=get(handles.EdgSV_popupmenu,'Value');
    end
end
EC.edg.size_ratio=str2double(get(handles.EdgSRR_edit,'String'));
EC.edg.size_abs=get(handles.EdgS_checkbox,'Value');
EC.edg.color_abs=get(handles.EdgC_checkbox,'Value');
if get(handles.EdgCS_radiobutton,'Value')==1
    EC.edg.color=1;
    EC.edg.CM(1,:)=get(handles.EdgCS_text,'BackgroundColor');
elseif get(handles.EdgCC_radiobutton,'Value')==1
    EC.edg.color=2;
    EC.edg.color_map=get(handles.EdgCC_popupmenu,'Value');
    switch EC.edg.color_map
        case 1
            EC.edg.CM=jet;
        case 2
            EC.edg.CM=hsv;
        case 3
            EC.edg.CM=hot;
        case 4
            EC.edg.CM=cool;
        case 5
            EC.edg.CM=spring;
        case 6
            EC.edg.CM=summer;
        case 7
            EC.edg.CM=autumn;
        case 8
            EC.edg.CM=winter;
        case 9
            EC.edg.CM=gray;
        case 10
            EC.edg.CM=bone;
        case 11
            EC.edg.CM=copper;
        case 12
            EC.edg.CM=pink;
        case 13
            EC.edg.CM=lines;
    end
elseif get(handles.EdgCT_radiobutton,'Value')==1
    EC.edg.color=3;
    EC.edg.color_threshold=str2double(get(handles.EdgCT_edit,'String'));
    EC.edg.CM(1,:)=get(handles.EdgCTH_text,'BackgroundColor');
    EC.edg.CM(64,:)=get(handles.EdgCTL_text,'BackgroundColor');
else
    EC.edg.color=4;
    EC.edg.color_distance=str2double(get(handles.EdgCD_edit,'String'));
    EC.edg.CM(1,:)=get(handles.EdgCDH_text,'BackgroundColor');
    EC.edg.CM(64,:)=get(handles.EdgCDL_text,'BackgroundColor');
end
EC.vol.display=get(handles.VolD_popupmenu,'Value');
EC.vol.pn=str2double(get(handles.VolPRn_edit,'String'));
EC.vol.px=str2double(get(handles.VolPRx_edit,'String'));
EC.vol.nn=str2double(get(handles.VolNRn_edit,'String'));
EC.vol.nx=str2double(get(handles.VolNRx_edit,'String'));
EC.vol.null=get(handles.VolNCS_text,'BackgroundColor');
EC.vol.color_map=get(handles.VolC_popupmenu,'Value');
EC.img.width=str2double(get(handles.ImgPW_edit,'String'));
EC.img.height=str2double(get(handles.ImgPH_edit,'String'));
EC.img.dpi=str2double(get(handles.ImgD_edit,'String'));

% --- Executes on button press in Cancel_button.
function Cancel_button_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global FLAG
FLAG.EC=0;
close(findobj('Tag','EP_fig'));

% --- Executes on slider movement.
function MshA_slider_Callback(hObject, eventdata, handles)
% hObject    handle to MshA_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = get(handles.MshA_slider, 'value');
set(handles.MshA_edit, 'string', num2str(val,'%6.2f'));
ChangeFlag(handles);


% --- Executes during object creation, after setting all properties.
function MshA_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MshA_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function MshA_edit_Callback(hObject, eventdata, handles)
% hObject    handle to MshA_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MshA_edit as text
%        str2double(get(hObject,'String')) returns contents of MshA_edit as a double
val=str2double(get(handles.MshA_edit,'String'));
set(handles.MshA_slider,'Value',val);
ChangeFlag(handles);


% --- Executes during object creation, after setting all properties.
function MshA_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MshA_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in NodCC_popupmenu.
function NodCC_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to NodCC_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns NodCC_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from NodCC_popupmenu
ChangeFlag(handles);


% --- Executes during object creation, after setting all properties.
function NodCC_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NodCC_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in NodCM_pushbutton.
function NodCM_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to NodCM_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
H=BrainNet_ModuleColor;
ChangeFlag(handles);



% --- Executes on slider movement.
function NodCT_slider_Callback(hObject, eventdata, handles)
% hObject    handle to NodCT_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = get(handles.NodCT_slider, 'value');
set(handles.NodCT_edit, 'string', num2str(val,'%6.2f'));
ChangeFlag(handles);


% --- Executes during object creation, after setting all properties.
function NodCT_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NodCT_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function NodCT_edit_Callback(hObject, eventdata, handles)
% hObject    handle to NodCT_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NodCT_edit as text
%        str2double(get(hObject,'String')) returns contents of NodCT_edit as a double
val=str2double(get(handles.NodCT_edit,'String'));
set(handles.NodCT_slider,'Value',val);
ChangeFlag(handles);


% --- Executes during object creation, after setting all properties.
function NodCT_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NodCT_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function NodST_slider_Callback(hObject, eventdata, handles)
% hObject    handle to NodST_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = get(handles.NodST_slider, 'value');
set(handles.NodST_edit, 'string', num2str(val,'%6.2f'));
ChangeFlag(handles);


% --- Executes during object creation, after setting all properties.
function NodST_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NodST_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function NodST_edit_Callback(hObject, eventdata, handles)
% hObject    handle to NodST_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NodST_edit as text
%        str2double(get(hObject,'String')) returns contents of NodST_edit as a double
val=str2double(get(handles.NodST_edit,'String'));
set(handles.NodST_slider,'Value',val);
ChangeFlag(handles);


% --- Executes during object creation, after setting all properties.
function NodST_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NodST_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function NodSVR_slider_Callback(hObject, eventdata, handles)
% hObject    handle to NodSVR_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = get(handles.NodSVR_slider, 'value');
set(handles.NodSVR_edit, 'string', num2str(val,'%6.2f'));
ChangeFlag(handles);


% --- Executes during object creation, after setting all properties.
function NodSVR_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NodSVR_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function NodSVR_edit_Callback(hObject, eventdata, handles)
% hObject    handle to NodSVR_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NodSVR_edit as text
%        str2double(get(hObject,'String')) returns contents of NodSVR_edit as a double
val=str2double(get(handles.NodSVR_edit,'String'));
set(handles.NodSVR_slider,'Value',val);
ChangeFlag(handles);


% --- Executes during object creation, after setting all properties.
function NodSVR_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NodSVR_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NodSS_edit_Callback(hObject, eventdata, handles)
% hObject    handle to NodSS_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NodSS_edit as text
%        str2double(get(hObject,'String')) returns contents of NodSS_edit as a double
ChangeFlag(handles);


% --- Executes during object creation, after setting all properties.
function NodSS_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NodSS_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in NodSV_popupmenu.
function NodSV_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to NodSV_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns NodSV_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from NodSV_popupmenu
global surf
if get(handles.NodSV_popupmenu,'Value')==2
    if min(surf.sphere(:,5))<1 || max(surf.sphere(:,5))>10
        msgbox('The size inputed may exceed the proper range, and will be adjusted!','Warning','warn');
    end
end
ChangeFlag(handles);


% --- Executes during object creation, after setting all properties.
function NodSV_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NodSV_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function NodDT_slider_Callback(hObject, eventdata, handles)
% hObject    handle to NodDT_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = get(handles.NodDT_slider, 'value');
set(handles.NodDT_edit, 'string', num2str(val,'%6.2f'));
ChangeFlag(handles);


% --- Executes during object creation, after setting all properties.
function NodDT_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NodDT_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function NodDT_edit_Callback(hObject, eventdata, handles)
% hObject    handle to NodDT_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NodDT_edit as text
%        str2double(get(hObject,'String')) returns contents of NodDT_edit as a double
val=str2double(get(handles.NodDT_edit,'String'));
set(handles.NodDT_slider,'Value',val);
ChangeFlag(handles);


% --- Executes during object creation, after setting all properties.
function NodDT_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NodDT_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in NodDT_popupmenu.
function NodDT_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to NodDT_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns NodDT_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from NodDT_popupmenu
global EC
global surf
str = get(hObject, 'string');
n = get(hObject, 'value');
switch str{n}
    case 'Size'
        if EC.nod.draw_threshold>max(surf.sphere(:,5)) || EC.nod.draw_threshold<min(surf.sphere(:,5))
            EC.nod.draw_threshold=min(surf.sphere(:,5));
        end
        set(handles.NodDT_slider,'Max',max(surf.sphere(:,5)));
        set(handles.NodDT_slider,'Min',min(surf.sphere(:,5)));
        set(handles.NodDT_slider,'Value',EC.nod.draw_threshold);
        set(handles.NodDT_edit,'String',num2str(EC.nod.draw_threshold,'%6.2f'));
    case 'Color'
        if EC.nod.draw_threshold>max(surf.sphere(:,4)) || EC.nod.draw_threshold<min(surf.sphere(:,4))
            EC.nod.draw_threshold=min(surf.sphere(:,4));
        end
        set(handles.NodDT_slider,'Max',max(surf.sphere(:,4)));
        set(handles.NodDT_slider,'Min',min(surf.sphere(:,4)));
        set(handles.NodDT_slider,'Value',EC.nod.draw_threshold);
        set(handles.NodDT_edit,'String',num2str(EC.nod.draw_threshold,'%6.2f'));
end
ChangeFlag(handles);



% --- Executes during object creation, after setting all properties.
function NodDT_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NodDT_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when Lbl_panel is resized.
function Lbl_panel_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to Lbl_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in LblA_radiobutton.
function LblA_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to LblA_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LblA_radiobutton


% --- Executes on button press in LblN_radiobutton.
function LblN_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to LblN_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LblN_radiobutton


% --- Executes on button press in LblT_radiobutton.
function LblT_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to LblT_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LblT_radiobutton


% --- Executes on slider movement.
function LblT_slider_Callback(hObject, eventdata, handles)
% hObject    handle to LblT_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = get(handles.LblT_slider, 'value');
set(handles.LblT_edit, 'string', num2str(val,'%6.2f'));
ChangeFlag(handles);


% --- Executes during object creation, after setting all properties.
function LblT_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LblT_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function LblT_edit_Callback(hObject, eventdata, handles)
% hObject    handle to LblT_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LblT_edit as text
%        str2double(get(hObject,'String')) returns contents of LblT_edit as a double
val=str2double(get(handles.LblT_edit,'String'));
set(handles.LblT_slider,'Value',val);
ChangeFlag(handles);



% --- Executes during object creation, after setting all properties.
function LblT_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LblT_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in LblT_popupmenu.
function LblT_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to LblT_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns LblT_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from LblT_popupmenu
global EC
global surf
str = get(hObject, 'string');
n = get(hObject, 'value');
switch str{n}
    case 'Size'
        if EC.lbl_threshold>max(surf.sphere(:,5)) || EC.lbl_threshold<min(surf.sphere(:,5))
            EC.lbl_threshold=min(surf.sphere(:,5));
        end
        set(handles.LblT_slider,'Max',max(surf.sphere(:,5)));
        set(handles.lblT_slider,'Min',min(surf.sphere(:,5)));
        set(handles.LblT_slider,'Value',EC.lbl_threshold);
        set(handles.LblT_edit,'String',num2str(EC.lbl_threshold,'%6.2f'));
    case 'Color'
        if EC.lbl_threshold>max(surf.sphere(:,4)) || EC.lbl_threshold<min(surf.sphere(:,4))
            EC.lbl_threshold=min(surf.sphere(:,4));
        end
        set(handles.LblT_slider,'Max',max(surf.sphere(:,4)));
        set(handles.lblT_slider,'Min',min(surf.sphere(:,4)));
        set(handles.LblT_slider,'Value',EC.lbl_threshold);
        set(handles.LblT_edit,'String',num2str(EC.lbl_threshold,'%6.2f'));
end
ChangeFlag(handles);



% --- Executes during object creation, after setting all properties.
function LblT_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LblT_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in EdgCC_popupmenu.
function EdgCC_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to EdgCC_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns EdgCC_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from EdgCC_popupmenu
ChangeFlag(handles);


% --- Executes during object creation, after setting all properties.
function EdgCC_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EdgCC_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function EdgCT_slider_Callback(hObject, eventdata, handles)
% hObject    handle to EdgCT_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = get(handles.EdgCT_slider, 'value');
set(handles.EdgCT_edit, 'string', num2str(val,'%6.2f'));
ChangeFlag(handles);


% --- Executes during object creation, after setting all properties.
function EdgCT_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EdgCT_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function EdgCT_edit_Callback(hObject, eventdata, handles)
% hObject    handle to EdgCT_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EdgCT_edit as text
%        str2double(get(hObject,'String')) returns contents of EdgCT_edit as a double
val=str2double(get(handles.EdgCT_edit,'String'));
set(handles.EdgCT_slider,'Value',val);
ChangeFlag(handles);


% --- Executes during object creation, after setting all properties.
function EdgCT_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EdgCT_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function EdgST_slider_Callback(hObject, eventdata, handles)
% hObject    handle to EdgST_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = get(handles.EdgST_slider, 'value');
set(handles.EdgST_edit, 'string', num2str(val,'%6.2f'));
ChangeFlag(handles);


% --- Executes during object creation, after setting all properties.
function EdgST_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EdgST_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function EdgST_edit_Callback(hObject, eventdata, handles)
% hObject    handle to EdgST_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EdgST_edit as text
%        str2double(get(hObject,'String')) returns contents of EdgST_edit as a double
val=str2double(get(handles.EdgST_edit,'String'));
set(handles.EdgST_slider,'Value',val);
ChangeFlag(handles);


% --- Executes during object creation, after setting all properties.
function EdgST_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EdgST_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EdgSRR_edit_Callback(hObject, eventdata, handles)
% hObject    handle to EdgSRR_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EdgSRR_edit as text
%        str2double(get(hObject,'String')) returns contents of EdgSRR_edit as a double
val=str2double(get(handles.EdgSRR_edit,'String'));
set(handles.EdgSRR_slider,'Value',val);
ChangeFlag(handles);


% --- Executes during object creation, after setting all properties.
function EdgSRR_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EdgSRR_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EdgSS_edit_Callback(hObject, eventdata, handles)
% hObject    handle to EdgSS_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EdgSS_edit as text
%        str2double(get(hObject,'String')) returns contents of EdgSS_edit as a double
ChangeFlag(handles);


% --- Executes during object creation, after setting all properties.
function EdgSS_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EdgSS_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in EdgSV_popupmenu.
function EdgSV_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to EdgSV_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns EdgSV_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from EdgSV_popupmenu
global surf
if get(handles.EdgSV_popupmenu,'Value')==2
    if min(min(surf.net))<0.2 || max(max(surf.net))>3
        msgbox('The size inputed may exceed the proper range, and will be adjusted!','Warning','warn');
    end
end
ChangeFlag(handles);


% --- Executes during object creation, after setting all properties.
function EdgSV_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EdgSV_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function EdgSRR_slider_Callback(hObject, eventdata, handles)
% hObject    handle to EdgSRR_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = get(handles.EdgSRR_slider, 'value');
set(handles.EdgSRR_edit, 'string', num2str(val,'%6.2f'));
ChangeFlag(handles);


% --- Executes during object creation, after setting all properties.
function EdgSRR_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EdgSRR_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function EdgDT_slider_Callback(hObject, eventdata, handles)
% hObject    handle to EdgDT_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global surf
val = get(handles.EdgDT_slider, 'value');
set(handles.EdgDT_edit, 'string', num2str(val,'%6.2f'));
switch get(handles.EdgDT_checkbox,'Value')
    case 0
        set(handles.EdgDS_slider,'Value',length(find(surf.net>val))/(size(surf.net,1)*size(surf.net,2)));
        set(handles.EdgDS_edit,'String',num2str(length(find(surf.net>val))/(size(surf.net,1)*size(surf.net,2)),'%6.2f'));
    case 1
        set(handles.EdgDS_slider,'Value',length(find(abs(surf.net)>val))/(size(surf.net,1)*size(surf.net,2)));
        set(handles.EdgDS_edit,'String',num2str(length(find(abs(surf.net)>val))/(size(surf.net,1)*size(surf.net,2)),'%6.2f'));
end
ChangeFlag(handles);



% --- Executes during object creation, after setting all properties.
function EdgDT_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EdgDT_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function EdgDT_edit_Callback(hObject, eventdata, handles)
% hObject    handle to EdgDT_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EdgDT_edit as text
%        str2double(get(hObject,'String')) returns contents of EdgDT_edit as a double
global surf
val=str2double(get(handles.EdgDT_edit,'String'));
set(handles.EdgDT_slider,'Value',val);
set(handles.EdgDS_slider,'Value',length(find(surf.net>EC.edg.draw_threshold))/(size(surf.net,1)*size(surf.net,2)));
set(handles.EdgDS_edit,'String',num2str(length(find(surf.net>EC.edg.draw_threshold))/(size(surf.net,1)*size(surf.net,2)),'%6.2f'));
ChangeFlag(handles);


% --- Executes during object creation, after setting all properties.
function EdgDT_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EdgDT_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in EdgDT_popupmenu.
function EdgDT_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to EdgDT_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns EdgDT_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from EdgDT_popupmenu


% --- Executes during object creation, after setting all properties.
function EdgDT_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EdgDT_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ImgPW_edit_Callback(hObject, eventdata, handles)
% hObject    handle to ImgPW_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ImgPW_edit as text
%        str2double(get(hObject,'String')) returns contents of ImgPW_edit as a double
global EC
switch get(handles.ImgDW_popupmenu,'Value')
    case 1
        set(handles.ImgDW_edit,'String',num2str(str2double(get(handles.ImgPW_edit,'String'))/str2double(get(handles.ImgD_edit,'String'))*2.54,'%6.2f'));
    case 2
        set(handles.ImgDW_edit,'String',num2str(str2double(get(handles.ImgPW_edit,'String'))/str2double(get(handles.ImgD_edit,'String')),'%6.2f'));
end
if get(handles.ImgC_checkbox,'Value')==1
    set(handles.ImgPH_edit,'String',num2str(str2double(get(handles.ImgPW_edit,'String'))*EC.img.height/EC.img.width,'%5d'));
    set(handles.ImgDH_edit,'String',num2str(str2double(get(handles.ImgDW_edit,'String'))*EC.img.height/EC.img.width,'%6.2f'));
end


% --- Executes during object creation, after setting all properties.
function ImgPW_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImgPW_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ImgPH_edit_Callback(hObject, eventdata, handles)
% hObject    handle to ImgPH_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ImgPH_edit as text
%        str2double(get(hObject,'String')) returns contents of ImgPH_edit as a double
global EC
switch get(handles.ImgDW_popupmenu,'Value')
    case 1
        set(handles.ImgDH_edit,'String',num2str(str2double(get(handles.ImgPH_edit,'String'))/str2double(get(handles.ImgD_edit,'String'))*2.54,'%6.2f'));
    case 2
        set(handles.ImgDH_edit,'String',num2str(str2double(get(handles.ImgPH_edit,'String'))/str2double(get(handles.ImgD_edit,'String')),'%6.2f'));
end
if get(handles.ImgC_checkbox,'Value')==1
    set(handles.ImgPW_edit,'String',num2str(str2double(get(handles.ImgPH_edit,'String'))*EC.img.width/EC.img.height,'%5d'));
    set(handles.ImgDW_edit,'String',num2str(str2double(get(handles.ImgDH_edit,'String'))*EC.img.width/EC.img.height,'%6.2f'));
end


% --- Executes during object creation, after setting all properties.
function ImgPH_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImgPH_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ImgD_edit_Callback(hObject, eventdata, handles)
% hObject    handle to ImgD_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ImgD_edit as text
%        str2double(get(hObject,'String')) returns contents of ImgD_edit as a double
switch get(handles.ImgDW_popupmenu,'Value')
    case 1
        set(handles.ImgPW_edit,'String',num2str(str2double(get(handles.ImgDW_edit,'String'))*str2double(get(handles.ImgD_edit,'String'))/2.54,'%5.0f'));
        set(handles.ImgPH_edit,'String',num2str(str2double(get(handles.ImgDH_edit,'String'))*str2double(get(handles.ImgD_edit,'String'))/2.54,'%5.0f'));
    case 2
        set(handles.ImgPW_edit,'String',num2str(str2double(get(handles.ImgDW_edit,'String'))*str2double(get(handles.ImgD_edit,'String')),'%5.0f'));
        set(handles.ImgPH_edit,'String',num2str(str2double(get(handles.ImgDH_edit,'String'))*str2double(get(handles.ImgD_edit,'String')),'%5.0f'));
end




% --- Executes during object creation, after setting all properties.
function ImgD_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImgD_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over BakC_text.
function BakC_text_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to BakC_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
c=uisetcolor('Select Color');
if length(c)==3
    set(handles.BakC_text,'BackgroundColor',c);
    ChangeFlag(handles);
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over MshC_text.
function MshC_text_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to MshC_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
c=uisetcolor('Select Color');
if length(c)==3
    set(handles.MshC_text,'BackgroundColor',c);
    ChangeFlag(handles);
end


% --- Executes when selected object is changed in NodD_panel.
function NodD_panel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in NodD_panel
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if get(handles.NodDA_radiobutton,'Value')==1
    set(handles.NodDT_slider,'Enable','off');
    set(handles.NodDT_edit,'Enable','off');
    set(handles.NodDT_popupmenu,'Enable','off');
else
    set(handles.NodDT_slider,'Enable','on');
    set(handles.NodDT_edit,'Enable','on');
    set(handles.NodDT_popupmenu,'Enable','on');
end
ChangeFlag(handles);


% --- Executes when selected object is changed in NodS_panel.
function NodS_panel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in NodS_panel
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if get(handles.NodSS_radiobutton,'Value')==1
    set(handles.NodSS_edit,'Enable','on');
    set(handles.NodSV_popupmenu,'Enable','off');
    set(handles.NodST_slider,'Enable','off');
    set(handles.NodST_edit,'Enable','off');
    set(handles.NodST_checkbox,'Enable','off');
elseif get(handles.NodSV_radiobutton,'Value')==1
    set(handles.NodSS_edit,'Enable','off');
    set(handles.NodSV_popupmenu,'Enable','on');
    set(handles.NodST_checkbox,'Enable','on');
    if get(handles.NodST_checkbox,'Value')==1
        set(handles.NodST_slider,'Enable','on');
        set(handles.NodST_edit,'Enable','on');
    else
        set(handles.NodST_slider,'Enable','off');
        set(handles.NodST_edit,'Enable','off');
    end
end
ChangeFlag(handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over NodCS_text.
function NodCS_text_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to NodCS_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
c=uisetcolor('Select Color');
if length(c)==3
    set(handles.NodCS_text,'BackgroundColor',c);
    ChangeFlag(handles);
end


% --- Executes when selected object is changed in NodC_panel.
function NodC_panel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in NodC_panel
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if get(handles.NodCS_radiobutton,'Value')==1
    set(handles.NodCC_popupmenu,'Enable','off');
    set(handles.NodCM_pushbutton,'Enable','off');
    set(handles.NodCT_slider,'Enable','off');
    set(handles.NodCT_edit,'Enable','off');
elseif get(handles.NodCC_radiobutton,'Value')==1
    set(handles.NodCC_popupmenu,'Enable','on');
    set(handles.NodCM_pushbutton,'Enable','off');
    set(handles.NodCT_slider,'Enable','off');
    set(handles.NodCT_edit,'Enable','off');
elseif get(handles.NodCM_radiobutton,'Value')==1
    set(handles.NodCC_popupmenu,'Enable','off');
    set(handles.NodCM_pushbutton,'Enable','on');
    set(handles.NodCT_slider,'Enable','off');
    set(handles.NodCT_edit,'Enable','off');
else
    set(handles.NodCC_popupmenu,'Enable','off');
    set(handles.NodCM_pushbutton,'Enable','off');
    set(handles.NodCT_slider,'Enable','on');
    set(handles.NodCT_edit,'Enable','on');
end
ChangeFlag(handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over NodCTH_text.
function NodCTH_text_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to NodCTH_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
c=uisetcolor('Select Color');
if length(c)==3
    set(handles.NodCTH_text,'BackgroundColor',c);
    ChangeFlag(handles);
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over NodCTL_text.
function NodCTL_text_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to NodCTL_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
c=uisetcolor('Select Color');
if length(c)==3
    set(handles.NodCTL_text,'BackgroundColor',c);
    ChangeFlag(handles);
end


% --- Executes when selected object is changed in Lbl_panel.
function Lbl_panel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in Lbl_panel
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if get(handles.LblA_radiobutton,'Value')==1 || get(handles.LblN_radiobutton,'Value')==1
    set(handles.LblT_slider,'Enable','off');
    set(handles.LblT_edit,'Enable','off');
    set(handles.LblT_popupmenu,'Enable','off');
else
    set(handles.LblT_slider,'Enable','on');
    set(handles.LblT_edit,'Enable','on');
    set(handles.LblT_popupmenu,'Enable','on');
end
ChangeFlag(handles);



% --- Executes when selected object is changed in EdgD_panel.
function EdgD_panel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in EdgD_panel
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if get(handles.EdgDA_radiobutton,'Value')==1
    set(handles.EdgDT_slider,'Enable','off');
    set(handles.EdgDT_edit,'Enable','off');
    set(handles.EdgDS_slider,'Enable','off');
    set(handles.EdgDS_edit,'Enable','off');
    set(handles.EdgDT_checkbox,'Enable','off');
else
    set(handles.EdgDT_slider,'Enable','on');
    set(handles.EdgDT_edit,'Enable','on');
    set(handles.EdgDS_slider,'Enable','on');
    set(handles.EdgDS_edit,'Enable','on');
    set(handles.EdgDT_checkbox,'Enable','on');
end
ChangeFlag(handles);



% --- Executes when selected object is changed in EdgS_uipanel.
function EdgS_uipanel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in EdgS_uipanel
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if get(handles.EdgSS_radiobutton,'Value')==1
    set(handles.EdgSS_edit,'Enable','on');
    set(handles.EdgSV_popupmenu,'Enable','off');
    set(handles.EdgST_slider,'Enable','off');
    set(handles.EdgST_edit,'Enable','off');
    set(handles.EdgST_checkbox,'Enable','off');
elseif get(handles.EdgSV_radiobutton,'Value')==1
    set(handles.EdgSS_edit,'Enable','off');
    set(handles.EdgSV_popupmenu,'Enable','on');
    set(handles.EdgST_checkbox,'Enable','on');
    if get(handles.EdgST_checkbox,'Value')==1
        set(handles.EdgST_slider,'Enable','on');
        set(handles.EdgST_edit,'Enable','on');
    else
        set(handles.EdgST_slider,'Enable','off');
        set(handles.EdgST_edit,'Enable','off');
    end
end
ChangeFlag(handles);


% --- Executes when selected object is changed in EdgC_panel.
function EdgC_panel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in EdgC_panel
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if get(handles.EdgCS_radiobutton,'Value')==1
    set(handles.EdgCC_popupmenu,'Enable','off');
    set(handles.EdgCT_slider,'Enable','off');
    set(handles.EdgCT_edit,'Enable','off');
    set(handles.EdgCD_slider,'Enable','off');
    set(handles.EdgCD_edit,'Enable','off');
elseif get(handles.EdgCC_radiobutton,'Value')==1
    set(handles.EdgCC_popupmenu,'Enable','on');
    set(handles.EdgCT_slider,'Enable','off');
    set(handles.EdgCT_edit,'Enable','off');
    set(handles.EdgCD_slider,'Enable','off');
    set(handles.EdgCD_edit,'Enable','off');
elseif get(handles.EdgCT_radiobutton,'Value')==1
    set(handles.EdgCC_popupmenu,'Enable','off');
    set(handles.EdgCT_slider,'Enable','on');
    set(handles.EdgCT_edit,'Enable','on');
    set(handles.EdgCD_slider,'Enable','off');
    set(handles.EdgCD_edit,'Enable','off');
else
    set(handles.EdgCC_popupmenu,'Enable','off');
    set(handles.EdgCT_slider,'Enable','off');
    set(handles.EdgCT_edit,'Enable','off');
    set(handles.EdgCD_slider,'Enable','on');
    set(handles.EdgCD_edit,'Enable','on');
end
ChangeFlag(handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over EdgCS_text.
function EdgCS_text_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to EdgCS_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
c=uisetcolor('Select Color');
if length(c)==3
    set(handles.EdgCS_text,'BackgroundColor',c);
end
ChangeFlag(handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over EdgCTH_text.
function EdgCTH_text_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to EdgCTH_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
c=uisetcolor('Select Color');
if length(c)==3
    set(handles.EdgCTH_text,'BackgroundColor',c);
    ChangeFlag(handles);
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over EdgCTL_text.
function EdgCTL_text_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to EdgCTL_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
c=uisetcolor('Select Color');
if length(c)==3
    set(handles.EdgCTL_text,'BackgroundColor',c);
    ChangeFlag(handles);
end


% --- Executes on button press in Load_button.
function Load_button_Callback(hObject, eventdata, handles)
% hObject    handle to Load_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname]=uigetfile({'*.mat','MAT-files (*.mat)'},'Load Configuration');
if isequal(filename,0) || isequal(pathname,0)
    return;
else
    fpath=fullfile(pathname,filename);
    load(fpath);
    Initialization(handles);
    msgbox('Option Loaded!','Success','help');
end

% --- Executes on button press in Save_button.
function Save_button_Callback(hObject, eventdata, handles)
% hObject    handle to Save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global EC
[filename,pathname]=uiputfile({'*.mat','MAT-files (*.mat)'},'Save Option');
if isequal(filename,0) || isequal(pathname,0)
    return;
else
    fpath=fullfile(pathname,filename);
    GetValue(handles);
    save(fpath,'EC');
    msgbox('Option Saved!','Success','help');
end


% --- Executes on button press in LblF_button.
function LblF_button_Callback(hObject, eventdata, handles)
% hObject    handle to LblF_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global EC
if isempty(EC.lbl_font)
    s=uisetfont;
else
    s=uisetfont(EC.lbl_font);
end
if isstruct(s)
    EC.lbl_font=s;
    ChangeFlag(handles);
end


% --- Executes on slider movement.
function EdgCD_slider_Callback(hObject, eventdata, handles)
% hObject    handle to EdgCD_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = get(handles.EdgCD_slider, 'value');
set(handles.EdgCD_edit, 'string', num2str(val,'%6.2f'));
ChangeFlag(handles);


% --- Executes during object creation, after setting all properties.
function EdgCD_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EdgCD_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function EdgCD_edit_Callback(hObject, eventdata, handles)
% hObject    handle to EdgCD_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EdgCD_edit as text
%        str2double(get(hObject,'String')) returns contents of EdgCD_edit as a double
val=str2double(get(handles.EdgCD_edit,'String'));
set(handles.EdgCD_slider,'Value',val);
ChangeFlag(handles);


% --- Executes during object creation, after setting all properties.
function EdgCD_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EdgCD_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over EdgCDH_text.
function EdgCDH_text_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to EdgCDH_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
c=uisetcolor('Select Color');
if length(c)==3
    set(handles.EdgCDH_text,'BackgroundColor',c);
    ChangeFlag(handles);
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over EdgCDL_text.
function EdgCDL_text_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to EdgCDL_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
c=uisetcolor('Select Color');
if length(c)==3
    set(handles.EdgCDL_text,'BackgroundColor',c);
    ChangeFlag(handles);
end


% --- Executes on slider movement.
function EdgDS_slider_Callback(hObject, eventdata, handles)
% hObject    handle to EdgDS_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global surf
val = get(handles.EdgDS_slider, 'value');
set(handles.EdgDS_edit, 'string', num2str(val,'%6.2f'));
sq=surf.nsph^2;
switch get(handles.EdgDT_checkbox,'Value')
    case 0
        temp=sort(reshape(surf.net,sq,1),'descend');
    case 1
        temp=sort(reshape(abs(surf.net),sq,1),'descend');
end
index=int32(val*sq);
if index<1
    index=1;
elseif index>sq
    index=sq;
end
set(handles.EdgDT_slider,'Value',temp(index));
set(handles.EdgDT_edit,'String',num2str(temp(index),'%6.2f'));
ChangeFlag(handles);


% --- Executes during object creation, after setting all properties.
function EdgDS_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EdgDS_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function EdgDS_edit_Callback(hObject, eventdata, handles)
% hObject    handle to EdgDS_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EdgDS_edit as text
%        str2double(get(hObject,'String')) returns contents of EdgDS_edit as a double
global surf
val = str2double(get(handles.EdgDS_edit, 'String'));
set(handles.EdgDS_slider, 'Value', val);
sq=surf.nsph^2;
switch get(handles.EdgDT_checkbox,'Value')
    case 0
        temp=sort(reshape(surf.net,sq,1),'descend');
    case 1
        temp=sort(reshape(abs(surf.net),sq,1),'descend');
end
index=int32(val*sq);
if index<1
    index=1;
elseif index>sq
    index=sq;
end
set(handles.EdgDT_slider,'Value',temp(index));
set(handles.EdgDT_edit,'String',num2str(temp(index),'%6.2f'));
ChangeFlag(handles);


% --- Executes during object creation, after setting all properties.
function EdgDS_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EdgDS_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on EdgDS_edit and none of its controls.
function EdgDS_edit_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to EdgDS_edit (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in EdgS_checkbox.
function EdgS_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to EdgS_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of EdgS_checkbox
global EC
global surf
switch get(handles.EdgS_checkbox,'Value')
    case 0
        if EC.edg.size_threshold>max(max(surf.net)) || EC.edg.size_threshold<min(min(surf.net))
            EC.edg.size_threshold=min(min(surf.net));
        end
        set(handles.EdgST_slider,'Max',max(max(surf.net)));
        set(handles.EdgST_slider,'Min',min(min(surf.net))-0.001);
        set(handles.EdgST_slider,'Value',EC.edg.size_threshold);
        set(handles.EdgST_edit,'String',num2str(EC.edg.size_threshold,'%6.2f'));
    case 1
        if EC.edg.size_threshold>max(max(abs(surf.net))) || EC.edg.size_threshold<min(min(abs(surf.net)))
            EC.edg.size_threshold=min(min(abs(surf.net)));
        end
        set(handles.EdgST_slider,'Max',max(max(abs(surf.net))));
        set(handles.EdgST_slider,'Min',min(min(abs(surf.net)))-0.001);
        set(handles.EdgST_slider,'Value',EC.edg.size_threshold);
        set(handles.EdgST_edit,'String',num2str(EC.edg.size_threshold,'%6.2f'));
end
ChangeFlag(handles);


% --- Executes on button press in EdgDT_checkbox.
function EdgDT_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to EdgDT_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of EdgDT_checkbox
global EC
global surf
switch get(handles.EdgDT_checkbox,'Value')
    case 0
        if EC.edg.draw_threshold>max(max(surf.net)) || EC.edg.draw_threshold<min(min(surf.net))
            EC.edg.draw_threshold=min(min(surf.net));
        end
        set(handles.EdgDT_slider,'Max',max(max(surf.net)));
        set(handles.EdgDT_slider,'Min',min(min(surf.net))-0.001);
        set(handles.EdgDT_slider,'Value',EC.edg.draw_threshold);
        set(handles.EdgDT_edit,'String',num2str(EC.edg.draw_threshold,'%6.2f'));
        set(handles.EdgDS_slider,'Value',length(find(surf.net>EC.edg.draw_threshold))/(size(surf.net,1)*size(surf.net,2)));
        set(handles.EdgDS_edit,'String',num2str(length(find(surf.net>EC.edg.draw_threshold))/(size(surf.net,1)*size(surf.net,2)),'%6.2f'));
    case 1
        if EC.edg.draw_threshold>max(max(abs(surf.net))) || EC.edg.draw_threshold<min(min(abs(surf.net)))
            EC.edg.draw_threshold=min(min(abs(surf.net)));
        end
        set(handles.EdgDT_slider,'Max',max(max(abs(surf.net))));
        set(handles.EdgDT_slider,'Min',min(min(abs(surf.net)))-0.001);
        set(handles.EdgDT_slider,'Value',EC.edg.draw_threshold);
        set(handles.EdgDT_edit,'String',num2str(EC.edg.draw_threshold,'%6.2f'));
        set(handles.EdgDS_slider,'Value',length(find(abs(surf.net)>EC.edg.draw_threshold))/(size(surf.net,1)*size(surf.net,2)));
        set(handles.EdgDS_edit,'String',num2str(length(find(abs(surf.net)>EC.edg.draw_threshold))/(size(surf.net,1)*size(surf.net,2)),'%6.2f'));
end
ChangeFlag(handles);


% --- Executes on button press in EdgC_checkbox.
function EdgC_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to EdgC_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of EdgC_checkbox
global EC
global surf
switch get(handles.EdgC_checkbox,'Value')
    case 0
        if EC.edg.color_threshold>max(max(surf.net)) || EC.edg.color_threshold<min(min(surf.net))
            EC.edg.color_threshold=min(min(surf.net));
        end
        set(handles.EdgCT_slider,'Max',max(max(surf.net)));
        set(handles.EdgCT_slider,'Min',min(min(surf.net))-0.001);
        set(handles.EdgCT_slider,'Value',EC.edg.color_threshold);
        set(handles.EdgCT_edit,'String',num2str(EC.edg.color_threshold,'%6.2f'));
    case 1
        if EC.edg.color_threshold>max(max(abs(surf.net))) || EC.edg.color_threshold<min(min(abs(surf.net)))
            EC.edg.color_threshold=min(min(abs(surf.net)));
        end
        set(handles.EdgCT_slider,'Max',max(max(abs(surf.net))));
        set(handles.EdgCT_slider,'Min',min(min(abs(surf.net)))-0.001);
        set(handles.EdgCT_slider,'Value',EC.edg.color_threshold);
        set(handles.EdgCT_edit,'String',num2str(EC.edg.color_threshold,'%6.2f'));
end
ChangeFlag(handles);


% --- Executes when selected object is changed in Lot_panel.
function Lot_panel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in Lot_panel
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if get(handles.LotP_radiobutton,'Value')==1
    set(handles.LotPS_radiobutton,'Enable','on');
    set(handles.LotPA_radiobutton,'Enable','on');
    set(handles.LotPC_radiobutton,'Enable','on');
else
    set(handles.LotPS_radiobutton,'Enable','off');
    set(handles.LotPA_radiobutton,'Enable','off');
    set(handles.LotPC_radiobutton,'Enable','off');
end
ChangeFlag(handles);

function ChangeFlag(handles)
global FLAG
FLAG.EC_change=1;
set(handles.Apply_button,'Enable','on');


% --- Executes on button press in Apply_button.
function Apply_button_Callback(hObject, eventdata, handles)
% hObject    handle to Apply_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global FLAG
FLAG.EC_change=0;
set(handles.Apply_button,'Enable','off');
GetValue(handles);
FLAG.EC=2;
uiresume(gcbf);




% --- Executes when selected object is changed in LotP_panel.
function LotP_panel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in LotP_panel
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
ChangeFlag(handles);


% --- Executes on button press in ImgC_checkbox.
function ImgC_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to ImgC_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ImgC_checkbox



function ImgDH_edit_Callback(hObject, eventdata, handles)
% hObject    handle to ImgDH_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ImgDH_edit as text
%        str2double(get(hObject,'String')) returns contents of ImgDH_edit as a double
global EC
switch get(handles.ImgDH_popupmenu,'Value')
    case 1
        set(handles.ImgPH_edit,'String',num2str(str2double(get(handles.ImgDH_edit,'String'))*str2double(get(handles.ImgD_edit,'String'))/2.54,'%5.0f'));
    case 2
        set(handles.ImgPH_edit,'String',num2str(str2double(get(handles.ImgDH_edit,'String'))*str2double(get(handles.ImgD_edit,'String')),'%5.0f'));
end
if get(handles.ImgC_checkbox,'Value')==1
    set(handles.ImgDW_edit,'String',num2str(str2double(get(handles.ImgDH_edit,'String'))*EC.img.width/EC.img.height,'%6.2f'));
    set(handles.ImgPW_edit,'String',num2str(str2double(get(handles.ImgPH_edit,'String'))*EC.img.width/EC.img.height,'%5.0f'));
end



% --- Executes during object creation, after setting all properties.
function ImgDH_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImgDH_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ImgDW_edit_Callback(hObject, eventdata, handles)
% hObject    handle to ImgDW_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ImgDW_edit as text
%        str2double(get(hObject,'String')) returns contents of ImgDW_edit as a double
global EC
switch get(handles.ImgDW_popupmenu,'Value')
    case 1
        set(handles.ImgPW_edit,'String',num2str(str2double(get(handles.ImgDW_edit,'String'))*str2double(get(handles.ImgD_edit,'String'))/2.54,'%5.0f'));
    case 2
        set(handles.ImgPW_edit,'String',num2str(str2double(get(handles.ImgDW_edit,'String'))*str2double(get(handles.ImgD_edit,'String')),'%5.0f'));
end
if get(handles.ImgC_checkbox,'Value')==1
    set(handles.ImgDH_edit,'String',num2str(str2double(get(handles.ImgDW_edit,'String'))*EC.img.height/EC.img.width,'%6.2f'));
    set(handles.ImgPH_edit,'String',num2str(str2double(get(handles.ImgPW_edit,'String'))*EC.img.height/EC.img.width,'%5.0f'));
end


% --- Executes during object creation, after setting all properties.
function ImgDW_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImgDW_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ImgDW_popupmenu.
function ImgDW_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to ImgDW_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ImgDW_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ImgDW_popupmenu
val=get(handles.ImgDW_popupmenu,'Value');
set(handles.ImgDH_popupmenu,'Value',val);
switch val
    case 1
        set(handles.ImgDW_edit,'String',num2str(str2double(get(handles.ImgPW_edit,'String'))/str2double(get(handles.ImgD_edit,'String'))*2.54,'%6.2f'));
        set(handles.ImgDH_edit,'String',num2str(str2double(get(handles.ImgPH_edit,'String'))/str2double(get(handles.ImgD_edit,'String'))*2.54,'%6.2f'));
    case 2
        set(handles.ImgDW_edit,'String',num2str(str2double(get(handles.ImgPW_edit,'String'))/str2double(get(handles.ImgD_edit,'String')),'%6.2f'));
        set(handles.ImgDH_edit,'String',num2str(str2double(get(handles.ImgPH_edit,'String'))/str2double(get(handles.ImgD_edit,'String')),'%6.2f'));
end



% --- Executes during object creation, after setting all properties.
function ImgDW_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImgDW_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ImgDH_popupmenu.
function ImgDH_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to ImgDH_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ImgDH_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ImgDH_popupmenu
val=get(handles.ImgDH_popupmenu,'Value');
set(handles.ImgDW_popupmenu,'Value',val);
switch val
    case 1
        set(handles.ImgDW_edit,'String',num2str(str2double(get(handles.ImgPW_edit,'String'))/str2double(get(handles.ImgD_edit,'String'))*2.54,'%6.2f'));
        set(handles.ImgDH_edit,'String',num2str(str2double(get(handles.ImgPH_edit,'String'))/str2double(get(handles.ImgD_edit,'String'))*2.54,'%6.2f'));
    case 2
        set(handles.ImgDW_edit,'String',num2str(str2double(get(handles.ImgPW_edit,'String'))/str2double(get(handles.ImgD_edit,'String')),'%6.2f'));
        set(handles.ImgDH_edit,'String',num2str(str2double(get(handles.ImgPH_edit,'String'))/str2double(get(handles.ImgD_edit,'String')),'%6.2f'));
end


% --- Executes during object creation, after setting all properties.
function ImgDH_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImgDH_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in NodST_checkbox.
function NodST_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to NodST_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of NodST_checkbox
if get(handles.NodST_checkbox,'Value')==1
    set(handles.NodST_slider,'Enable','on');
    set(handles.NodST_edit,'Enable','on');
else
    set(handles.NodST_slider,'Enable','off');
    set(handles.NodST_edit,'Enable','off');
end
ChangeFlag(handles);


% --- Executes on button press in EdgST_checkbox.
function EdgST_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to EdgST_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of EdgST_checkbox
if get(handles.EdgST_checkbox,'Value')==1
    set(handles.EdgST_slider,'Enable','on');
    set(handles.EdgST_edit,'Enable','on');
else
    set(handles.EdgST_slider,'Enable','off');
    set(handles.EdgST_edit,'Enable','off');
end
ChangeFlag(handles);


% --- Executes on selection change in VolC_popupmenu.
function VolC_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to VolC_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns VolC_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from VolC_popupmenu
ChangeFlag(handles);

% --- Executes during object creation, after setting all properties.
function VolC_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VolC_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function VolNRn_edit_Callback(hObject, eventdata, handles)
% hObject    handle to VolNRn_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VolNRn_edit as text
%        str2double(get(hObject,'String')) returns contents of VolNRn_edit as a double
ChangeFlag(handles);

% --- Executes during object creation, after setting all properties.
function VolNRn_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VolNRn_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function VolNRx_edit_Callback(hObject, eventdata, handles)
% hObject    handle to VolNRx_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VolNRx_edit as text
%        str2double(get(hObject,'String')) returns contents of VolNRx_edit as a double
ChangeFlag(handles);

% --- Executes during object creation, after setting all properties.
function VolNRx_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VolNRx_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function VolPRn_edit_Callback(hObject, eventdata, handles)
% hObject    handle to VolPRn_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VolPRn_edit as text
%        str2double(get(hObject,'String')) returns contents of VolPRn_edit as a double
ChangeFlag(handles);

% --- Executes during object creation, after setting all properties.
function VolPRn_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VolPRn_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function VolPRx_edit_Callback(hObject, eventdata, handles)
% hObject    handle to VolPRx_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VolPRx_edit as text
%        str2double(get(hObject,'String')) returns contents of VolPRx_edit as a double
ChangeFlag(handles);

% --- Executes during object creation, after setting all properties.
function VolPRx_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VolPRx_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in VolD_popupmenu.
function VolD_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to VolD_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns VolD_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from VolD_popupmenu
val=get(hObject,'Value');
switch val
    case 1
        set(handles.VolNR_text,'Enable','on');
        set(handles.VolNRn_edit,'Enable','on');
        set(handles.VolNRx_edit,'Enable','on');
        set(handles.VolPR_text,'Enable','on');
        set(handles.VolPRn_edit,'Enable','on');
        set(handles.VolPRx_edit,'Enable','on');
    case 2
        set(handles.VolNR_text,'Enable','off');
        set(handles.VolNRn_edit,'Enable','off');
        set(handles.VolNRx_edit,'Enable','off');
        set(handles.VolPR_text,'Enable','on');
        set(handles.VolPRn_edit,'Enable','on');
        set(handles.VolPRx_edit,'Enable','on');
    case 3
        set(handles.VolPR_text,'Enable','off');
        set(handles.VolPRn_edit,'Enable','off');
        set(handles.VolPRx_edit,'Enable','off');
        set(handles.VolNR_text,'Enable','on');
        set(handles.VolNRn_edit,'Enable','on');
        set(handles.VolNRx_edit,'Enable','on');
        set(handles.VolPR_text,'Enable','on');
end
ChangeFlag(handles);


% --- Executes during object creation, after setting all properties.
function VolD_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VolD_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over VolNCS_text.
function VolNCS_text_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to VolNCS_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
c=uisetcolor('Select Color');
if length(c)==3
    set(hObject,'BackgroundColor',c);
    ChangeFlag(handles);
end
