function varargout = insertText(varargin)
% function varargout = insertText(varargin)
% INSERTTEXT M-file for insertText.fig
% INSERTTEXT, by itself, creates a new INSERTTEXT or raises the existing
%      singleton*.
%
%      H = INSERTTEXT returns the handle to a new INSERTTEXT or the handle to
%      the existing singleton*.
%
%      INSERTTEXT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INSERTTEXT.M with the given input arguments.
%
%      INSERTTEXT('Property','Value',...) creates a new INSERTTEXT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before insertText_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to insertText_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help insertText

% Last Modified by GUIDE v2.5 21-Sep-2007 20:29:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @insertText_OpeningFcn, ...
                   'gui_OutputFcn',  @insertText_OutputFcn, ...
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


% --- Executes just before insertText is made visible.
function insertText_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to insertText (see VARARGIN)

% Choose default command line output for insertText
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes insertText wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = insertText_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox.
function listbox_Callback(hObject, eventdata, handles)
% hObject    handle to listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox


% --- Executes during object creation, after setting all properties.
function listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonInsert.
function pushbuttonInsert_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonInsert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dictory = get(handles.editDir,'String');
rekursiv = get(handles.checkboxSubDirs,'Value');
textFile = get(handles.editText,'String');

fid = fopen([pwd,filesep,textFile],'r');
tline = fgetl(fid);
curLine = 1;
while ischar(tline)
    text{curLine} = tline;
    curLine = curLine+1;
    tline = fgetl(fid);
end
            

processDictory(dictory,text,rekursiv==1);


function processDictory(dictory,text,rekursiv)

mfiles = what(dictory);
if ~isempty(mfiles)
    mfiles = mfiles.m;
    for curMfile = 1:size(mfiles)
        fprintf([mfiles{curMfile},'\n']);
        if ~isempty(text)
            % read the hole file first
            fid = fopen([pwd,filesep,dictory,filesep,mfiles{curMfile}],'r+');
            tline = fgetl(fid);
            curLine = 1;
            clear fileData
            while ischar(tline)
                fileData{curLine} = tline;
                curLine = curLine+1;
                tline = fgetl(fid);
            end
            curLine = 1;
            while (curLine<=size(fileData,2)) && ~isempty(fileData{curLine}) && ...
                    (strcmp(fileData{curLine}(1),'%') || ~isempty(findstr(fileData{curLine},'function')))
                curLine = curLine+1;
            end
            inserLine = curLine-1;
            % write back
            frewind(fid);
            if inserLine == 0
                for j = 1 : size(text,2)
                    fprintf(fid,'%s\n',text{j});
                end
            end
            for curLine = 1 : size(fileData,2)
                fprintf(fid,'%s\n',fileData{curLine});
                if curLine == inserLine
                    for j = 1 : size(text,2)-1
                        fprintf(fid,'%s\n',text{j});
                    end
                    fprintf(fid,'%s',text{end});
                end
            end
            fclose(fid);
        end   
    end
end

if rekursiv
    dictories = dir(dictory);
    if ~isempty(dictories)
        for curDir = 1:size(dictories)
            if dictories(curDir).isdir && isempty(findstr('.',dictories(curDir).name))
                dirName = [dictory,filesep,dictories(curDir).name];
                strDirName = dirName;
                strDirName(findstr(strDirName,'\')) = '/';
                fprintf(['\n',strDirName,'\n ___________________________\n']);
                processDictory(dirName,text,rekursiv);
            end
        end
    end
end






function editDir_Callback(hObject, eventdata, handles)
% hObject    handle to editDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDir as text
%        str2double(get(hObject,'String')) returns contents of editDir as a double


% --- Executes during object creation, after setting all properties.
function editDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxSubDirs.
function checkboxSubDirs_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxSubDirs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxSubDirs





function editText_Callback(hObject, eventdata, handles)
% hObject    handle to editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editText as text
%        str2double(get(hObject,'String')) returns contents of editText as a double


% --- Executes during object creation, after setting all properties.
function editText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


