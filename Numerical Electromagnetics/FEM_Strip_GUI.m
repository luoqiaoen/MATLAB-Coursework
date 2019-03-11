function varargout = FEM_Strip_GUI(varargin)
% FEM_STRIP_GUI M-file for FEM_Strip_GUI.fig
%      FEM_STRIP_GUI, by itself, creates a new FEM_STRIP_GUI or raises the existing
%      singleton*.
%
%      H = FEM_STRIP_GUI returns the handle to a new FEM_STRIP_GUI or the handle to
%      the existing singleton*.
%
%      FEM_STRIP_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FEM_STRIP_GUI.M with the given input arguments.
%
%      FEM_STRIP_GUI('Property','Value',...) creates a new FEM_STRIP_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FEM_Strip_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FEM_Strip_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FEM_Strip_GUI

% Last Modified by GUIDE v2.5 08-Apr-2016 11:58:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FEM_Strip_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @FEM_Strip_GUI_OutputFcn, ...
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


% --- Executes just before FEM_Strip_GUI is made visible.
function FEM_Strip_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FEM_Strip_GUI (see VARARGIN)

% Choose default command line output for FEM_Strip_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FEM_Strip_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FEM_Strip_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function d_y_Callback(hObject, eventdata, handles)
% hObject    handle to d_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of d_y as text
%        str2double(get(hObject,'String')) returns contents of d_y as a double


% --- Executes during object creation, after setting all properties.
function d_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to d_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function d_x_Callback(hObject, eventdata, handles)
% hObject    handle to d_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of d_x as text
%        str2double(get(hObject,'String')) returns contents of d_x as a double


% --- Executes during object creation, after setting all properties.
function d_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to d_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function strip_v_Callback(hObject, eventdata, handles)
% hObject    handle to strip_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of strip_v as text
%        str2double(get(hObject,'String')) returns contents of strip_v as a double


% --- Executes during object creation, after setting all properties.
function strip_v_CreateFcn(hObject, eventdata, handles)
% hObject    handle to strip_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function s_length_Callback(hObject, eventdata, handles)
% hObject    handle to s_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of s_length as text
%        str2double(get(hObject,'String')) returns contents of s_length as a double


% --- Executes during object creation, after setting all properties.
function s_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to s_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in simulate.
function simulate_Callback(hObject, eventdata, handles)
% hObject    handle to simulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

d_x =str2double(get(handles.d_x,'string'));
d_y=str2double(get(handles.d_y,'string'));
s_length=str2double(get(handles.s_length,'string'));
s_vol=str2double(get(handles.strip_v,'string'));


[phi, Ex, Ey] = FEM_Strip(d_x,d_y,s_length,s_vol);

axes(handles.plot)
hold on;
contourf(phi,0:s_vol/50:s_vol);
quiver(Ex,Ey);
colorbar;
axis image;
hold off;

% --- Executes on button press in quit.
function quit_Callback(hObject, eventdata, handles)
% hObject    handle to quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf)
