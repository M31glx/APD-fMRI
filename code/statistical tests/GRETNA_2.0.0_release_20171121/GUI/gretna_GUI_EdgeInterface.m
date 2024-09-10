function varargout = gretna_GUI_EdgeInterface(varargin)
% GRETNA_GUI_EDGEINTERFACE MATLAB code for gretna_GUI_EdgeInterface.fig
%      GRETNA_GUI_EDGEINTERFACE, by itself, creates a new GRETNA_GUI_EDGEINTERFACE or raises the existing
%      singleton*.
%
%      H = GRETNA_GUI_EDGEINTERFACE returns the handle to a new GRETNA_GUI_EDGEINTERFACE or the handle to
%      the existing singleton*.
%
%      GRETNA_GUI_EDGEINTERFACE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GRETNA_GUI_EDGEINTERFACE.M with the given input arguments.
%
%      GRETNA_GUI_EDGEINTERFACE('Property','Value',...) creates a new GRETNA_GUI_EDGEINTERFACE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gretna_GUI_EdgeInterface_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gretna_GUI_EdgeInterface_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gretna_GUI_EdgeInterface

% Last Modified by GUIDE v2.5 23-Feb-2017 17:41:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gretna_GUI_EdgeInterface_OpeningFcn, ...
                   'gui_OutputFcn',  @gretna_GUI_EdgeInterface_OutputFcn, ...
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


% --- Executes just before gretna_GUI_EdgeInterface is made visible.
function gretna_GUI_EdgeInterface_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gretna_GUI_EdgeInterface (see VARARGIN)
handles.Group1Cells={};
handles.Group2Cells={};
handles.CovCells={};
handles.CurDir=pwd;
set(handles.OutputDirEntry, 'String', pwd);
% Choose default command line output for gretna_GUI_EdgeInterface
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gretna_GUI_EdgeInterface wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gretna_GUI_EdgeInterface_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Group2Button.
function Group2Button_Callback(hObject, eventdata, handles)
% hObject    handle to Group2Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[Name, Path]=uigetfile({'*.txt;*.mat','Brain Network Matrix (*.txt;*.mat)';'*.*', 'All Files (*.*)';}, ...
            'Pick brain network matrix' , 'MultiSelect' , 'On', handles.CurDir);
if isnumeric(Name)
    return
end
handles.CurDir=Path;

if ischar(Name)
    Name={Name};
end
Name=Name';

NameCell=cellfun(@(name) fullfile(Path, name), Name, 'UniformOutput', false);
handles.Group2Cells=[handles.Group2Cells; NameCell];
handles=GenDataListbox(handles);
guidata(hObject, handles);


% --- Executes on button press in Group2Button.
function Group1Button_Callback(hObject, eventdata, handles)
% hObject    handle to Group2Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[Name, Path]=uigetfile({'*.txt;*.mat','Brain Network Matrix (*.txt;*.mat)';'*.*', 'All Files (*.*)';}, ...
            'Pick brain network matrix' , 'MultiSelect' , 'On', handles.CurDir);
if isnumeric(Name)
    return
end
handles.CurDir=Path;

if ischar(Name)
    Name={Name};
end
Name=Name';

NameCell=cellfun(@(name) fullfile(Path, name), Name, 'UniformOutput', false);
handles.Group1Cells=[handles.Group1Cells; NameCell];
handles=GenDataListbox(handles);
guidata(hObject, handles);

% --- Executes on selection change in GroupListbox.
function GroupListbox_Callback(hObject, eventdata, handles)
% hObject    handle to GroupListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(gcf , 'SelectionType') , 'normal')
    return;
elseif strcmp(get(gcf , 'SelectionType') , 'open')
    SelectedValue=get(handles.GroupListbox , 'Value');
    if SelectedValue==0
        return
    end
    DataText=get(handles.GroupListbox , 'String'); 
    
    Info=DataText{SelectedValue};
    if exist(Info, 'file')==2
        Index=0;
        for i=1:numel(handles.Group1Cells)
            if strcmpi(handles.Group1Cells{i}, Info)
                Index=i;
                break
            end
        end
        if Index
            handles.Group1Cells(Index)=[];
        end
        Index=0;
        for i=1:numel(handles.Group2Cells)
            if strcmpi(handles.Group1Cells{i}, Info)
                Index=i;
                break
            end
        end
        if Index
            handles.Group1Cells(Index)=[];
        end
        handles=GenDataListbox(handles);
        guidata(hObject, handles);
        return
    end
    Tok=regexp(Info, 'Group\d: (.*?)_.*', 'tokens');
    
    if isempty(Tok)
        return
    end
    FFlag=Tok{1}{1};
    
    Match=regexp(DataText, 'Group\d: (.*?)_.*');
    for i=SelectedValue:-1:1
        if isempty(Match{i})
            break
        end
    end
    PathName=DataText{i};
    
    if strcmpi(FFlag, 'MAT')
        TempMat=load(PathName);
        Tok=regexp(Info, 'Group\d: MAT_.*_VAR_(.*):.*', 'tokens');
        NAME=Tok{1}{1};
        if ~isfield(TempMat, NAME)
            Tok=regexp(Info, 'Group\d: MAT_.*_VAR_(.*)_CELL_.*:.*', 'tokens');
            if isempty(Tok)
                Tok=regexp(Info, 'Group\d: MAT_.*_VAR_(.*)_STRUCT_.*:.*', 'tokens');
                NAME=Tok{1}{1};
            else
                NAME=Tok{1}{1};
            end
        end
        VAR=TempMat.(NAME);
        if isnumeric(VAR)
            Net=VAR;
            figure('Name', sprintf('%s(VAR: %s)', PathName, NAME),...
                'Numbertitle', 'Off'),
            imagesc(Net);
        elseif iscell(VAR)
            Tok=regexp(Info, 'Group\d: MAT_.*_VAR_.*_CELL_(.*):.*', 'tokens');
            NUM=str2num(Tok{1}{1});
            Net=VAR{NUM};
            figure('Name', sprintf('%s(VAR: %s --> %.4d)', PathName, NAME, NUM),...
                'Numbertitle', 'Off'),
            imagesc(Net);
        elseif isstruct(VAR)
            Tok=regexp(Info, 'Group\d: MAT_.*_VAR_.*_STRUCT_(.*):.*', 'tokens');
            SUB=Tok{1}{1};
            Net=VAR.(SUB);
            figure('Name', sprintf('%s(VAR: %s --> %s)', PathName, NAME, SUB),...
                'Numbertitle', 'Off'),
            imagesc(Net);
        end
    elseif strcmpi(FFlag, 'TXT')
        Net=load(PathName);
        figure('Name', sprintf('%s', PathName), 'Numbertitle' , 'Off'),
        imagesc(Net);
    end
end
% Hints: contents = cellstr(get(hObject,'String')) returns GroupListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GroupListbox


% --- Executes during object creation, after setting all properties.
function GroupListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GroupListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in TypePopup.
function TypePopup_Callback(hObject, eventdata, handles)
% hObject    handle to TypePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Value=get(handles.TypePopup, 'Value');
switch Value
    case 1 %Average
        CorrectLabelString='Threshold Type';
        CorrectPopupString={'Similarity threshold';'Sparsity';'LANS'};
        CorrectPupupValue=2;
        TLFlag='Off';
        TPFlag='Off';
        PLFlag='Off';
        PEFlag='Off';
        PLabelString='Threshold';
        PEntryString='0.05';
        G2Flag='Off';
        CFlag='Off';
        BFlag='Off';
        handles.Group2Cells={};
    case 2 %Backbone
        CorrectLabelString='Threshold Type';
        CorrectPopupString={'Edge Probability'};
        CorrectPupupValue=1;
        TLFlag='On';
        TPFlag='On';
        PLFlag='On';
        PEFlag='On';        
        PLabelString='Threshold';
        PEntryString='0.25';
        G2Flag='Off';
        CFlag='Off';
        BFlag='Off';
        handles.Group2Cells={};            
    case 3 %One Sample T-test
        CorrectLabelString='Correct Method';
        CorrectPopupString={'None';'FDR';'Bonferroni';'NBS'};
        CorrectPupupValue=1;
        TLFlag='On';
        TPFlag='On';
        PLFlag='On';
        PEFlag='On';         
        PLabelString='Edge P';
        PEntryString='0.05';
        G2Flag='Off'; 
        CFlag='On';
        BFlag='On';
        handles.Group2Cells={};
    case 4 %Two Sample T-test
        CorrectLabelString='Correct Method';
        CorrectPopupString={'None';'FDR';'Bonferroni';'NBS'};
        CorrectPupupValue=1;
        TLFlag='On';
        TPFlag='On';
        PLFlag='On';
        PEFlag='On';         
        PLabelString='Edge P';
        PEntryString='0.05';
        G2Flag='On';
        CFlag='On';
        BFlag='Off';
        handles.Group2Cells={};            
end
set(handles.BaseLabel, 'Visible', BFlag);
set(handles.BaseEntry, 'Visible', BFlag);
set(handles.Group2Button, 'Visible', G2Flag);
set(handles.CovAddButton, 'Enable', CFlag);
set(handles.CovRemoveButton, 'Enable', CFlag);
set(handles.CovListbox, 'Enable', CFlag, 'String', '');
handles.CovCells={};

handles.Group1Cells={};
handles.Group2Cells={};
set(handles.GroupListbox, 'String', '');
%handles=GenDataListbox(handles);

set(handles.CorrectLabel, 'String', CorrectLabelString, 'Visible', TLFlag);
set(handles.CorrectPopup, 'String', CorrectPopupString, 'Visible', TPFlag,...
    'Value', CorrectPupupValue);
set(handles.PLabel, 'String', PLabelString, 'Visible', PLFlag);
set(handles.PEntry, 'String', PEntryString, 'Visible', PEFlag);
set(handles.IterLabel,  'Visible', 'Off');
set(handles.IterEntry,  'Visible', 'Off');
set(handles.CPLabel,    'Visible', 'Off');
set(handles.CPEntry,    'Visible', 'Off');
set(handles.NMskLabel,  'Visible', 'Off');
set(handles.NMskEntry,  'Visible', 'Off', 'String', '');
set(handles.NMskButton, 'Visible', 'Off');
guidata(hObject, handles);

% Hints: contents = cellstr(get(hObject,'String')) returns TypePopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TypePopup


% --- Executes during object creation, after setting all properties.
function TypePopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TypePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PEntry_Callback(hObject, eventdata, handles)
% hObject    handle to PEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PEntry as text
%        str2double(get(hObject,'String')) returns contents of PEntry as a double


% --- Executes during object creation, after setting all properties.
function PEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in CorrectPopup.
function CorrectPopup_Callback(hObject, eventdata, handles)
% hObject    handle to CorrectPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Value=get(handles.CorrectPopup, 'Value');
TypeValue=get(handles.TypePopup, 'Value');
if TypeValue==1
    switch Value
        case 1 %Similarity threshold
            PLabelString='Threshold';
            Flag='Off';
        case 2 %Sparity
            PLabelString='Threshold';
            Flag='Off';
        case 3 %LANS
            PLabelString='Alpha';
            Flag='Off';
    end
elseif TypeValue==2
    switch Value
        case 1 % Edge Probability
            PLabelString='Threshold';
            Flag='Off';
        case 2 % Sign Test
            PLabelString='Edge P';
            Flag='Off'; 
    end
else
    switch Value
        case 1 %None
            PLabelString='Edge P';
            Flag='Off';
        case 2 %FDR
            PLabelString='FDR Q';
            Flag='Off';
        case 3 %Bonferroni
            PLabelString='Edge P';
            Flag='Off';
        case 4 %NBS
            PLabelString='Edge P';
            Flag='On';
    end
end
set(handles.PLabel,     'String',  PLabelString);
set(handles.IterLabel,  'VIsible', Flag);
set(handles.IterEntry,  'Visible', Flag);
set(handles.CPLabel,    'Visible', Flag);
set(handles.CPEntry,    'Visible', Flag);
set(handles.NMskLabel,  'Visible', Flag);
set(handles.NMskEntry,  'Visible', Flag);
set(handles.NMskButton, 'Visible', Flag);
% Hints: contents = cellstr(get(hObject,'String')) returns CorrectPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CorrectPopup


% --- Executes during object creation, after setting all properties.
function CorrectPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CorrectPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function OutputDirEntry_Callback(hObject, eventdata, handles)
% hObject    handle to OutputDirEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of OutputDirEntry as text
%        str2double(get(hObject,'String')) returns contents of OutputDirEntry as a double


% --- Executes during object creation, after setting all properties.
function OutputDirEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OutputDirEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in OutputDirButton.
function OutputDirButton_Callback(hObject, eventdata, handles)
% hObject    handle to OutputDirButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Path=uigetdir(handles.CurDir, 'Pick Output Directory');
if isnumeric(Path)
    return
end
handles.CurDir=fileparts(Path);
set(handles.OutputDirEntry, 'String', Path);
guidata(hObject, handles);



function PrefixEntry_Callback(hObject, eventdata, handles)
% hObject    handle to PrefixEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PrefixEntry as text
%        str2double(get(hObject,'String')) returns contents of PrefixEntry as a double


% --- Executes during object creation, after setting all properties.
function PrefixEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PrefixEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ComputeButton.
function ComputeButton_Callback(hObject, eventdata, handles)
% hObject    handle to ComputeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.Group1Cells)
    errordlg('Please import a group of network matrix');
    return;
end

fprintf('Computing...\n');

Type=get(handles.TypePopup, 'Value');
if Type==4
    if isempty(handles.Group2Cells)
        errordlg('Please import another group of network matrix');
        return;
    end
end

[MatrixGroup1, AliasList1]=GetGroupData(handles.Group1Cells);
fprintf('Group1:\n');
for i=1:numel(AliasList1)
    fprintf('\t%s\n', AliasList1{i});
end
    
OutputDir=get(handles.OutputDirEntry, 'String');
if isempty(OutputDir)
    OutputDir=fileparts(handles.CurDir);
end
Prefix=get(handles.PrefixEntry, 'String');
switch Type
    case 1 %Average
        AllMatrix=zeros([size(MatrixGroup1{1}),numel(MatrixGroup1)]);
        for i=1:numel(MatrixGroup1)
            Matrix=MatrixGroup1{i};
            Matrix=Matrix-diag(diag(Matrix));
            AllMatrix(:,:,i)=Matrix;
        end
        Matrix=mean(AllMatrix, 3);

        save(fullfile(OutputDir, [Prefix, '_Avg.txt']), 'Matrix', '-ASCII', '-DOUBLE', '-TABS');        
        fprintf('\nFinished.\n');
    case 2 % Structural Backbone
        AllRawMatrix=zeros([size(MatrixGroup1{1}),numel(MatrixGroup1)]);
        AllNumMatrix=zeros(size(AllRawMatrix));
        for i=1:numel(MatrixGroup1)
            Matrix=MatrixGroup1{i};
            Matrix=Matrix-diag(diag(Matrix));
            AllRawMatrix(:,:,i)=Matrix;
            AllNumMatrix(:, :, i)=logical(Matrix);
        end
        NumMatrix=sum(AllNumMatrix, 3);
        SumMatrix=sum(AllRawMatrix, 3);
        
        ProbMatrix=NumMatrix./numel(MatrixGroup1);
        MeanMatrix=SumMatrix./NumMatrix;
        MeanMatrix(isnan(MeanMatrix))=0;
        % Edge Probability
        ThresType=get(handles.CorrectPopup, 'Value');
        if ThresType==1
            ThresValue=str2double(get(handles.PEntry, 'String'));
        end
        
        BackboneMap=(ProbMatrix>ThresValue).*MeanMatrix;
        
        save(fullfile(OutputDir, [Prefix, '_Backbone.txt']), 'BackboneMap', '-ASCII', '-DOUBLE', '-TABS');
        save(fullfile(OutputDir, [Prefix, '_Probability.txt']), 'ProbMatrix', '-ASCII', '-DOUBLE', '-TABS'); 
        
        fprintf('\nFinished.\n');
    case 3 % One Sample T-test
        AllMatrix=zeros([size(MatrixGroup1{1}),numel(MatrixGroup1)]);
        for i=1:numel(MatrixGroup1)
            Matrix=MatrixGroup1{i};
            Matrix=Matrix-diag(diag(Matrix));
            AllMatrix(:,:,i)=Matrix;
        end
        [n1, n2, n3]=size(AllMatrix);
        AllMatrix=reshape(AllMatrix, [], n3);
        MIndex=triu(true(n1, n2), 1);
        AllMatrix=AllMatrix(MIndex(:), :);
        
        Group1Matrix={AllMatrix'};
        CovCells=handles.CovCells;
        for i=1:numel(CovCells)
            CovCells{i}=load(CovCells{i});
        end
        Base=str2double(get(handles.BaseEntry, 'String'));
        [T, P]=gretna_TTest1(Group1Matrix, CovCells, Base);

        CorrectValue=get(handles.CorrectPopup, 'Value');
        PorQ=str2double(get(handles.PEntry, 'String'));
        if CorrectValue==1 %None
            PThrd=PorQ;
            Index=find(P<PThrd);
            if isempty(Index)
                msgbox('No Edge Left');
                fprintf('\n\tApplying Edge P: Done.\n');
                return
            end
            TThrd=min(abs(T(Index)));
            fprintf('\n\tApplying Edge P: Done.\n');                        
        elseif CorrectValue==2 %FDR
            PThrd=gretna_FDR(P, PorQ);
            if isempty(PThrd)
                msgbox('No Edge Left'); 
                fprintf('\n\tFDR: Done.\n');
                return
            end
            TThrd=min(abs(T(P<PThrd)));            
            fprintf('\n\tFDR: Done.\n');
        elseif CorrectValue==3 %Bonferroni
            N=length(P);
            PThrd=PorQ/N;
            Index=find(P<PThrd);
            if isempty(Index)
                msgbox('No Edge Left\n'); 
                fprintf('\n\Bonferroni: Done.\n');
                return
            end
            TThrd=min(abs(T(Index)));            
            fprintf('\n\tBonferroni: Done.\n');
        elseif CorrectValue==4 %NBS
            PThrd=PorQ;
            Index=find(P<PThrd);
            if isempty(Index)
                msgbox('No Edge Left'); 
                fprintf('\n\Applying Edge P: Done.\n');
                return
            end            
            TThrd=min(abs(T(Index)));
            CPThrd=str2num(get(handles.CPEntry, 'String'));
            
            % NBS
            NMskName=get(handles.NMskEntry, 'String');
            if isempty(NMskName)
                NMsk=true(n1, n2);
            else
                NMsk=load(NMskName);
            end
            TMap=zeros(n1*n2, 1);
            TMap(MIndex(:))=T;
            TMap=reshape(TMap, [n1, n2]);
            TMap=TMap+TMap';            
            
            PMap=zeros(n1*n2, 1);
            PMap(MIndex(:))=P;
            PMap=reshape(PMap, [n1, n2]);
            PMap=PMap+PMap';
            M=str2double(get(handles.IterEntry, 'String'));
            
            [Comnet, Comnet_P]=gretna_TTest1_NBS(Group1Matrix,...
                NMsk, CovCells, Base, PThrd, CPThrd, TMap, PMap, M);
            if isempty(Comnet)
                msgbox(sprintf('No Comnet Left under Edge P<%f, Comnet P < %f', PThrd, CPThrd));
            else
                save(fullfile(OutputDir, [Prefix, '_ComnetMat.mat']), 'Comnet', 'Comnet_P');
                save(fullfile(OutputDir, [Prefix, '_ComnetP.txt']),  'Comnet_P', '-ASCII',  '-DOUBLE', '-TABS');
                for c=1:size(Comnet, 1)
                    comnet=Comnet{c, 2};
                    save(fullfile(OutputDir, sprintf('%s_%s_Comnet_%s.txt', Prefix, Comnet{c, 1}, Comnet{c, 3})),  'comnet', '-ASCII',  '-DOUBLE', '-TABS');
                end
            end
            fprintf('\n\tNBS: Done.\n');
        end
        TMap=zeros(n1*n2, 1);
        PMap=zeros(n1*n2, 1);
        TMap(MIndex(:))=T;
        PMap(MIndex(:))=P;
        TMap=reshape(TMap, [n1, n2]);
        PMap=reshape(PMap, [n1, n2]);
        TMap=TMap+TMap';
        PMap=PMap+PMap';
                    
        save(fullfile(OutputDir, [Prefix, '_TNet.txt']),  'TMap', '-ASCII', '-DOUBLE', '-TABS');
        save(fullfile(OutputDir, [Prefix, '_PNet.txt']),  'PMap', '-ASCII', '-DOUBLE', '-TABS');        
        save(fullfile(OutputDir, [Prefix, '_TThrd.txt']), 'TThrd', '-ASCII', '-DOUBLE', '-TABS');
        save(fullfile(OutputDir, [Prefix, '_PThrd.txt']), 'PThrd', '-ASCII', '-DOUBLE', '-TABS');        
        fprintf('\nFinished.\n');
    case 4 %Two Sample T-test
        [MatrixGroup2, AliasList2]=GetGroupData(handles.Group2Cells);
        fprintf('Group2:\n');
        for i=1:numel(AliasList2)
            fprintf('\t%s\n', AliasList2{i});
        end        
        AllMatrix1=zeros([size(MatrixGroup1{1}),numel(MatrixGroup1)]);
        AllMatrix2=zeros([size(MatrixGroup2{1}),numel(MatrixGroup2)]);
        %Group1
        for i=1:numel(MatrixGroup1)
            Matrix=MatrixGroup1{i};
            Matrix=Matrix-diag(diag(Matrix));
            AllMatrix1(:,:,i)=Matrix;
        end
        %Group2
        for i=1:numel(MatrixGroup2)
            Matrix=MatrixGroup2{i};
            Matrix=Matrix-diag(diag(Matrix));
            AllMatrix2(:,:,i)=Matrix;
        end        
        [n11, n12, n13]=size(AllMatrix1);
        [n21, n22, n23]=size(AllMatrix2);
        %Group1
        AllMatrix1=reshape(AllMatrix1, [], n13);
        MIndex=triu(true(n11, n12), 1); 
        AllMatrix1=AllMatrix1(MIndex(:), :);
        %Group2
        AllMatrix2=reshape(AllMatrix2, [], n23);
        MIndex=triu(true(n21, n22), 1); 
        AllMatrix2=AllMatrix2(MIndex(:), :);
        
        GroupMatrix=cell(2, 1);
        GroupMatrix{1, 1}=AllMatrix1';
        GroupMatrix{2, 1}=AllMatrix2';
        CovCells=handles.CovCells;
        for i=1:numel(CovCells)
            CovCells{i}=load(CovCells{i});
        end
        [T, P]=gretna_TTest2(GroupMatrix, CovCells);

        CorrectValue=get(handles.CorrectPopup, 'Value');
        PorQ=str2double(get(handles.PEntry, 'String'));
        if CorrectValue==1 %None
            PThrd=PorQ;
            Index=find(P<PThrd);
            if isempty(Index)
                msgbox('No Edge Left');
                fprintf('\n\tApplying Edge P: Done.\n');
                return
            end
            TThrd=min(abs(T(Index)));
            fprintf('\n\tApplying Edge P: Done.\n');            
        elseif CorrectValue==2 %FDR
            PThrd=gretna_FDR(P, PorQ);
            if isempty(PThrd)
                msgbox('No Edge Left'); 
                fprintf('\n\tFDR: Done.\n');
                return
            end
            TThrd=min(abs(T(P<PThrd)));            
            fprintf('\n\tFDR: Done.\n');
        elseif CorrectValue==3 %Bonferroni
            N=length(P);
            PThrd=PorQ/N;
            Index=find(P<PThrd);
            if isempty(Index)
                msgbox('No Edge Left\n'); 
                fprintf('\n\Bonferroni: Done.\n');
                return
            end
            TThrd=min(abs(T(Index)));            
            fprintf('\n\tBonferroni: Done.\n');
        elseif CorrectValue==4 %NBS
            PThrd=PorQ;
            Index=find(P<PThrd);
            if isempty(Index)
                msgbox('No Edge Left'); 
                fprintf('\n\tApplying Edge P: Done.\n');
                return
            end            
            TThrd=min(abs(T(Index)));
            CPThrd=str2num(get(handles.CPEntry, 'String'));
            
            % NBS
            NMskName=get(handles.NMskEntry, 'String');
            if isempty(NMskName)
                NMsk=true(n11, n12);
            else
                NMsk=load(NMskName);
            end
            TMap=zeros(n11*n12, 1);
            TMap(MIndex(:))=T;
            TMap=reshape(TMap, [n11, n12]);
            TMap=TMap+TMap';            
            
            PMap=zeros(n11*n12, 1);
            PMap(MIndex(:))=P;
            PMap=reshape(PMap, [n11, n12]);
            PMap=PMap+PMap';
            M=str2double(get(handles.IterEntry, 'String'));
            
            [Comnet, Comnet_P]=gretna_TTest2_NBS(GroupMatrix,...
                NMsk, CovCells, PThrd, CPThrd, TMap, PMap, M);
            if isempty(Comnet)
                msgbox(sprintf('No Comnet Left under Edge P < %g, Comnet P < %g', PThrd, CPThrd));
            else
                save(fullfile(OutputDir, [Prefix, '_ComnetMat.mat']), 'Comnet', 'Comnet_P');
                save(fullfile(OutputDir, [Prefix, '_ComnetP.txt']),  'Comnet_P', '-ASCII',  '-DOUBLE', '-TABS');
                for c=1:size(Comnet, 1)
                    comnet=Comnet{c, 2};
                    save(fullfile(OutputDir, sprintf('%s_%s_Comnet_%s.txt', Prefix, Comnet{c, 1}, Comnet{c, 3})),  'comnet', '-ASCII',  '-DOUBLE', '-TABS');
                end
            end
            fprintf('\n\tNBS: Done.\n'); 
        end
        TMap=zeros(n11*n12, 1);
        PMap=zeros(n11*n12, 1);
        TMap(MIndex(:))=T;
        PMap(MIndex(:))=P;
        TMap=reshape(TMap, [n11, n12]);
        PMap=reshape(PMap, [n11, n12]);
        TMap=TMap+TMap';
        PMap=PMap+PMap';
   
        save(fullfile(OutputDir, [Prefix, '_TNet.txt']),  'TMap', '-ASCII',  '-DOUBLE', '-TABS');
        save(fullfile(OutputDir, [Prefix, '_PNet.txt']),  'PMap', '-ASCII',  '-DOUBLE', '-TABS');        
        save(fullfile(OutputDir, [Prefix, '_TThrd.txt']), 'TThrd', '-ASCII', '-DOUBLE', '-TABS');
        save(fullfile(OutputDir, [Prefix, '_PThrd.txt']), 'PThrd', '-ASCII', '-DOUBLE', '-TABS'); 
        fprintf('\nFinished.\n');
end

function [MatrixGroup, AliasList]=GetGroupData(GroupCells) 
MatrixGroup=[];
AliasList=[];
for i=1:size(GroupCells)
    [Path , Name , Ext]=fileparts(GroupCells{i});
    if strcmp(Ext , '.txt')
        TempMat=load([Path , filesep , Name , Ext]);
        MatrixGroup=[MatrixGroup ;...
            {TempMat}];
        AliasList=[AliasList;{sprintf('TXT_%s', Name)}];
    elseif strcmp(Ext , '.mat')
        TempStruct=load([Path , filesep , Name , Ext]);
        FieldNames=fieldnames(TempStruct);
        for j=1:numel(FieldNames)
            TempMat=TempStruct.(FieldNames{j});
            if iscell(TempMat)
                MatrixGroup=[MatrixGroup; TempMat];
                for k=1:numel(TempMat)
                    AliasList=[AliasList;{sprintf('MAT_%s_VAR_%s_CELL_%.4d',...
                        Name, FieldNames{j}, k)}];
                end
            elseif isstruct(TempMat)
                SubField=fieldnames(TempMat);
                for k=1:numel(SubField)
                    MatrixGroup=[MatrixGroup; {TempMat.(SubField{k})}];
                    AliasList=[AliasList;{sprintf('MAT_%s_VAR_%s_STRUCT_%s',...
                        Name, FieldNames{j}, SubField{k})}];                   
                end
            else
                MatrixGroup=[MatrixGroup;{TempMat}];
                AliasList=[AliasList;{sprintf('MAT_%s_VAR_%s', Name, FieldNames{j})}];
            end
        end
    end
end

function handles=GenDataListbox(handles)
    present_list=[];
    Value=0;
    if ~isempty(handles.Group1Cells)
        Value=1;
        for i=1:numel(handles.Group1Cells)
        	[Path , Name , Ext]=...
            	fileparts(handles.Group1Cells{i});
            if strcmp(Ext , '.txt')
            	TempMat=load([Path , filesep , Name , Ext]);
                present_list=[present_list ;...
                	{fullfile(Path , [Name , Ext])}];
                present_list=[present_list ; ...
                	{sprintf('Group1: TXT_%s: (%d -- %d)', Name,...
                    size(TempMat , 1), size(TempMat , 2))}];
            elseif strcmp(Ext , '.mat')
            	TempStruct=load([Path , filesep , Name , Ext]);
                FieldNames=fieldnames(TempStruct);
                present_list=[present_list;...
                    {fullfile(Path, [Name, Ext])}];

                for j=1:numel(FieldNames)
                    TempMat=TempStruct.(FieldNames{j});
                    if iscell(TempMat)
                        for k=1:numel(TempMat)
                            present_list=[present_list ; ...
                                {sprintf('Group1: MAT_%s_VAR_%s_CELL_%.4d: (%d -- %d)',...
                                Name, FieldNames{j}, k,...
                                size(TempMat{k}, 1), size(TempMat{k}, 2))}];
                        end
                    elseif isstruct(TempMat)
                        SubField=fieldnames(TempMat);
                        for k=1:numel(SubField)
                            present_list=[present_list ; ...
                                {sprintf('Group1: MAT_%s_VAR_%s_STRUCT_%s: (%d -- %d)',...
                                Name, FieldNames{j}, SubField{k},...
                                size(TempMat.(SubField{k}), 1), size(TempMat.(SubField{k}), 2))}];
                        end
                    else
                        present_list=[present_list ; ...
                            {sprintf('Group1: MAT_%s_VAR_%s: (%d -- %d)',...
                            Name, FieldNames{j},...
                            size(TempMat , 1), size(TempMat , 2))}];
                    end
                end
            end
        end
    end
    
    if ~isempty(handles.Group2Cells)
        Value=1;
        for i=1:numel(handles.Group2Cells)
        	[Path , Name , Ext]=...
            	fileparts(handles.Group2Cells{i});
            if strcmp(Ext , '.txt')
            	TempMat=load([Path , filesep , Name , Ext]);
                present_list=[present_list ;...
                	{fullfile(Path , [Name , Ext])}];
                present_list=[present_list ; ...
                	{sprintf('Group2: TXT_%s: (%d -- %d)', Name,...
                    size(TempMat , 1), size(TempMat , 2))}];
            elseif strcmp(Ext , '.mat')
            	TempStruct=load([Path , filesep , Name , Ext]);
                FieldNames=fieldnames(TempStruct);
                present_list=[present_list;...
                    {fullfile(Path, [Name, Ext])}];

                for j=1:numel(FieldNames)
                    TempMat=TempStruct.(FieldNames{j});
                    if iscell(TempMat)
                        for k=1:numel(TempMat)
                            present_list=[present_list ; ...
                                {sprintf('Group2: MAT_%s_VAR_%s_CELL_%.4d: (%d -- %d)',...
                                Name, FieldNames{j}, k,...
                                size(TempMat{k}, 1), size(TempMat{k}, 2))}];
                        end
                    elseif isstruct(TempMat)
                        SubField=fieldnames(TempMat);
                        for k=1:numel(SubField)
                            present_list=[present_list ; ...
                                {sprintf('Group2: MAT_%s_VAR_%s_STRUCT_%s: (%d -- %d)',...
                                Name, FieldNames{j}, SubField{k},...
                                size(TempMat.(SubField{k}), 1), size(TempMat.(SubField{k}), 2))}];
                        end
                    else
                        present_list=[present_list ; ...
                            {sprintf('Group2: MAT_%s_VAR_%s: (%d -- %d)',...
                            Name, FieldNames{j},...
                            size(TempMat , 1), size(TempMat , 2))}];
                    end
                end
            end
        end
    end
    
    set(handles.GroupListbox , 'String' , present_list, 'Value', Value);
    drawnow;


% --- Executes on selection change in NetCutPopup.
function NetCutPopup_Callback(hObject, eventdata, handles)
% hObject    handle to NetCutPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns NetCutPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from NetCutPopup


% --- Executes during object creation, after setting all properties.
function NetCutPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NetCutPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function IterEntry_Callback(hObject, eventdata, handles)
% hObject    handle to IterEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IterEntry as text
%        str2double(get(hObject,'String')) returns contents of IterEntry as a double


% --- Executes during object creation, after setting all properties.
function IterEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IterEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CovAddButton.
function CovAddButton_Callback(hObject, eventdata, handles)
% hObject    handle to CovAddButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[Name, Path]=uigetfile({'*.txt;*.csv;*.tsv','Text Covariates File (*.txt;*.csv;*.tsv)';'*.*', 'All Files (*.*)';},...
    'Pick the Text Covariates', 'MultiSelect','on');
if isnumeric(Name)
    return
end

if ischar(Name)
    Name={Name};
end
Name=Name';

PathCell=cellfun(@(name) fullfile(Path, name), Name, 'UniformOutput', false);

handles.CovCells=[handles.CovCells;PathCell];

AddString(handles.CovListbox, PathCell);
guidata(hObject, handles);

% --- Executes on button press in CovRemoveButton.
function CovRemoveButton_Callback(hObject, eventdata, handles)
% hObject    handle to CovRemoveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Value=get(handles.CovListbox, 'Value');
if Value==0
    return
end
handles.CovCells(Value)=[];
RemoveString(handles.CovListbox, Value);
guidata(hObject, handles);

function AddString(ListboxHandle, NewCell)
StringCell=get(ListboxHandle, 'String');
StringCell=[StringCell; NewCell];
set(ListboxHandle, 'String', StringCell, 'Value', numel(StringCell));

function RemoveString(ListboxHandle, Value)
StringCell=get(ListboxHandle, 'String');
StringCell(Value)=[];
if isempty(StringCell)
    Value=0;
end
if Value > numel(StringCell)
    Value=Value-1;
end
set(ListboxHandle, 'String', StringCell, 'Value', Value);

% --- Executes on selection change in CovListbox.
function CovListbox_Callback(hObject, eventdata, handles)
% hObject    handle to CovListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CovListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CovListbox


% --- Executes during object creation, after setting all properties.
function CovListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CovListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BaseEntry_Callback(hObject, eventdata, handles)
% hObject    handle to BaseEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BaseEntry as text
%        str2double(get(hObject,'String')) returns contents of BaseEntry as a double


% --- Executes during object creation, after setting all properties.
function BaseEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BaseEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CPEntry_Callback(hObject, eventdata, handles)
% hObject    handle to CPEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CPEntry as text
%        str2double(get(hObject,'String')) returns contents of CPEntry as a double


% --- Executes during object creation, after setting all properties.
function CPEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CPEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NMskEntry_Callback(hObject, eventdata, handles)
% hObject    handle to NMskEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NMskEntry as text
%        str2double(get(hObject,'String')) returns contents of NMskEntry as a double


% --- Executes during object creation, after setting all properties.
function NMskEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NMskEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in NMskButton.
function NMskButton_Callback(hObject, eventdata, handles)
% hObject    handle to NMskButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[Name, Path]=uigetfile({'*.txt;*.csv;*.tsv','Text Covariates File (*.txt;*.csv;*.tsv)';'*.*', 'All Files (*.*)';},...
    'Pick the Text Covariates');
if isnumeric(Name)
    return
end

FullPath=fullfile(Path, Name);
set(handles.NMskEntry, 'String', FullPath);
