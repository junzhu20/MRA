function varargout = MRA(varargin)
% MRA MATLAB code for MRA.fig
% Wirtten by Yiwei Hou as part of the paper 
% "Multi-resolution analysis enables fidelity-ensured computational super-resolution and denoising for fluorescence microscopy"
% Authority: Yiwei Hou, Peng Xi
% College of fzuture technology, Peking University
% For any questions, please contact: houyiwei@stu.pku.edu.cn; xipeng@pku.edu.cn
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MRA_OpeningFcn, ...
                   'gui_OutputFcn',  @MRA_OutputFcn, ...
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


% --- Executes just before MRA is made visible.
function MRA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MRA (see VARARGIN)

% Choose default command line output for MRA
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
global filename;
global pathname;
global bit;
global f_stack_old;

% UIWAIT makes MRA wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MRA_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Input.
function Input_Callback(hObject, eventdata, handles)
% hObject    handle to Input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
addpath('./3D_Framelet_Denoising')
addpath('./2D_MRA_Deconvolution')
addpath('./3D_Comlex Dual Tree_Denoising')
addpath('./Basic_Operation')
global filename;
global pathname;
global bit;
global f_stack_old;
[filename,pathname]=uigetfile({'*.*'},"Select image");
if isequal(filename,0)||isequal(pathname,0)
   errordlg("No image selected","Error");
else
pat='.tif';
datatype=contains(filename,pat);
datatype=double(datatype);
if datatype~=1
       errordlg("Please select .tif images","Error");
end
bit=16;
f_stack_old=readMTiffn([pathname,filename],bit);
I=f_stack_old(:,:,1);
axes(handles.axes1);
imshow(double(I),[]);
end


% --- Executes on button press in Process.
function Process_Callback(hObject, eventdata, handles)
% hObject    handle to Process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global filename;
global pathname;
global bit;
global f_stack_old;
f_stack=f_stack_old;
flag1=get(handles.Activation3D,'value');
flag2=get(handles.Activation2D,'value');
[d1,d2,d3]=size(f_stack);
[~,sparsity_curvelet] = Sparsitycal_fun(f_stack(:,:,1));

%% Error check
% no option selected
if (flag1+flag2)==0
    msgbox('Please select at least one option','Error')
    error('Please select at least one option')
end
% 2D data input while activate spatiotemporal denoising
if d3==1&&flag1==1
 msgbox('Only one frame input,spatiotemporal continuity denoising turned off','Warning')
 flag1=0;
end
% check PSF parameter
if isempty(str2num(get(handles.Wavelength,'str')))||isempty(str2num(get(handles.Pixelsize,'str')))||isempty(str2num(get(handles.NA,'str')))||isempty(str2num(get(handles.Scale,'str')))==1
    msgbox('Please check the parameter input in PSF parameter column, which should be numbers','Error')
    error('Please check the parameter input in PSF parameter column, which should be numbers')
end
% check MRA parameter
if isempty(str2num(get(handles.maxframe,'str')))||isempty(str2num(get(handles.lambda1,'str')))||isempty(str2num(get(handles.lambda_framelet,'str')))||isempty(str2num(get(handles.lambda_curvelet,'str')))==1
    msgbox('Please check the parameter input in 2D MRA deconvolution or spatiotemporal continuity denoising column, which should be numbers','Error')
    error('Please check the parameter input in 2D MRA deconvolution or spatiotemporal continuity denoising column, which should be numbers')
end
% Check if spatiotemporal denoising is activated before
pat='Spatiotemporal denoised_';
flag11=contains(filename,pat);
flag11=double(flag11);

%% Padding
% Padding reflecsive
num1=20;
m=d1;
n=d2; 
f0_stack=zeros(m+2*num1,n+2*num1,d3,'single');
for i=1:d3
f0_stack(num1+1:end-num1,num1+1:end-num1,i)=f_stack(:,:,i);
f1=flip(f_stack(:,:,i));
f0_stack(1:num1,num1+1:end-num1,i)=f1(end-num1+1:end,1:end);
f0_stack(end-num1+1:end,num1+1:end-num1,i)=f1(1:num1,1:end);
f2=flip(f_stack(:,:,i),2);
f0_stack(num1+1:end-num1,1:num1,i)=f2(1:end,end-num1+1:end);
f0_stack(num1+1:end-num1,end-num1+1:end,i)=f2(1:end,1:num1);
f0_stack(1:num1,1:num1,i)=flip(flip(f_stack(1:num1,1:num1,i),2)');
f0_stack(end-num1+1:end,end-num1+1:end,i)=flip(flip(f_stack(end-num1+1:end,end-num1+1:end,i),2)');
f0_stack(1:num1,end-num1+1:end,i)=flip(flip(f_stack(1:num1,end-num1+1:end,i),2)');
f0_stack(end-num1+1:end,1:num1,i)=flip(flip(f_stack(end-num1+1:end,1:num1,i),2)');
end
f1_stack=f0_stack;
%Padding zeros
[d1,d2,d3]=size(f1_stack);
num2=10;
f0_stack=zeros(d1+2*num2,d2+2*num2,d3,'single');
f0_stack(num2+1:end-num2,num2+1:end-num2,:)=f1_stack;
[d1,d2,d3]=size(f0_stack);
num=num1+num2;

%% Generate PSF
if flag2==1
    [d1,d2,d3]=size(f0_stack);
lambda=get(handles.Wavelength,'str');
pixelsize=get(handles.Pixelsize,'str');
NA=get(handles.NA,'str');
factor=get(handles.Scale,'str');
lambda=str2num(lambda);
pixelsize=str2num(pixelsize);
NA=str2num(NA);
factor=str2num(factor);
PSF = PSF_Generator(lambda,pixelsize,NA,min([d1,d2]),factor);
end


%% 3D Framelet Pre-denoising & 3D Complex Dual-tree Pre-denoising
tic
if flag1==1
maxframe=get(handles.maxframe,'str');
maxframe=str2num(maxframe);
sub_n=ceil(d3/maxframe);
lambda1=get(handles.lambda1,'str');
lambda1=str2num(lambda1);
lamda_bank=lambda1*ones(1,ceil(d3/maxframe));
lamda_bank=single(lamda_bank);

if d3>1
    h = waitbar(0,'Spatiotemporal continuity denoising processing');
if d3>maxframe
   for iter=1:sub_n
     if iter~=sub_n&&iter~=1
         f_stack_single=f0_stack(:,:,1+(iter-1)*maxframe-2:iter*maxframe+2);    
     elseif iter==sub_n
         f_stack_single=f0_stack(:,:,1+(iter-1)*maxframe-2:end);     
     else
         f_stack_single=f0_stack(:,:,1:maxframe+2);     
     end
    lambda1=lamda_bank(iter);
    lambda2=single(lambda1/3);
    gpu=0;
    [f_stack_single] = FrameletDenoising3D(f_stack_single,lambda1,gpu);
        s=sprintf('Spatiotemporal continuity denoising processing:%d',ceil((iter-0.5)/sub_n*100));
        waitbar((iter-0.5)/sub_n,h,[s '%']);
    [f_stack_single] = DualTreeDenoising3D(f_stack_single,lambda2);
        s=sprintf('Spatiotemporal continuity denoising processing:%d',ceil((iter)/sub_n*100));
        waitbar((iter)/sub_n,h,[s '%']);
    if iter==1
    f0_stack(:,:,1:maxframe)=f_stack_single(:,:,1:maxframe);
    elseif iter==sub_n
    f0_stack(:,:,1+(iter-1)*maxframe:end)=f_stack_single(:,:,3:end);
    else
    f0_stack(:,:,1+(iter-1)*maxframe:iter*maxframe)=f_stack_single(:,:,3:end-2);
    end
   end
else
    lambda2=single(lambda1/3);
    gpu=0;
    [f0_stack] = FrameletDenoising3D(f0_stack,lambda1,gpu);
 s=sprintf('Spatiotemporal continuity denoising processing:%d',ceil(0.5*100));
       waitbar(0.5,h,[s '%']);
    [f0_stack] = DualTreeDenoising3D(f0_stack,lambda2);
end
end
        figure(99)
        imshow(f0_stack(:,:,1),[])
        title('Frame1: Spatiotemporal denoised'); pause(0.5)
 s=sprintf('Spatiotemporal continuity denoising processing:%d',ceil(1*100));
       waitbar(1,h,[s '%']);
        close(h)
end
%% 2D MRA Deconvolution
if flag2==1
lamda_framelet=get(handles.lambda_framelet,'str');
lamda_curvelet=get(handles.lambda_curvelet,'str');
RLtime=get(handles.RLtimes,'str');
lamda_framelet=str2num(lamda_framelet);
lamda_curvelet=str2num(lamda_curvelet);
RLtime=str2num(RLtime);
fprintf('**2D MRA deconvolution starting**\n');
[f0_stack,fun_all] = MRADeconvolution2D(f0_stack,PSF,lamda_framelet,lamda_curvelet,RLtime,flag1,flag11);
toc
fprintf('**2D MRA deconvolution completed**\n');
end
% f_stack=f0_stack;
% Regather
f_stack=f0_stack(num+1:end-num,num+1:end-num,:);

%% Write file
outputname=[pathname,'MRA Results\'];
mkdir(outputname)
[outputname] = GenerateTimeFolder(outputname);
filepara=outputname;
if flag1==1&&flag2~=1
outputname=[outputname,['Spatiotemporal denoised_',filename]];
    filepara=[outputname,'Spatiotemporal denoised_'];
    paraflag=1;
elseif flag1~=1&&flag2==1
    outputname=[outputname,['2D MRA deconvolved_',filename]];
    filepara=[outputname,'2D MRA deconvolved_'];
    paraflag=2;
else
     outputname=[outputname,['Full MRA deconvolved_',filename]];
     filepara=[outputname,'Full MRA deconvolved_'];
     paraflag=3;
end
writeMTiffn(f_stack,outputname,bit)
axes(handles.axes2);
imshow(f_stack(:,:,1),[]);
%% Write parameter
fp=fopen([filepara,'.txt'],'a');
if flag1==1
    fprintf(fp,'%3s ', '****Spatiotemporal continuity denoising parameter****');    
    fprintf(fp, '\n');
    fprintf(fp,'%3s ', '3D Framelet thresholding=');
    fprintf(fp, '%f\t',lambda1);
    fprintf(fp, '\n');
    fprintf(fp,'%3s ', 'Maximum consecutive frame=');
    fprintf(fp, '%f\t',maxframe);
    fprintf(fp, '\n');
end
if flag2==1
    fprintf(fp,'%3s ', '****PSF parameter****');    
    fprintf(fp, '\n');
    fprintf(fp,'%3s ', 'Wavelength=');
    fprintf(fp, '%f\t', lambda);
    fprintf(fp, '\n');
    fprintf(fp,'%3s ', 'Pixelsize=');
    fprintf(fp, '%f\t',pixelsize);
    fprintf(fp, '\n');
    fprintf(fp,'%3s ', 'NA=');
    fprintf(fp, '%f\t',NA);
    fprintf(fp, '\n');
    fprintf(fp,'%3s ', 'factor=');
    fprintf(fp, '%f\t',factor);
    fprintf(fp, '\n');
    
    fprintf(fp,'%3s ', '****2D MRA deconvolution parameter****');    
    fprintf(fp, '\n');
    fprintf(fp,'%3s ', 'Framelet thresholding=');
    fprintf(fp, '%f\t',lamda_framelet);
    fprintf(fp, '\n');
    fprintf(fp,'%3s ', 'Curvelet thresholding=');
    fprintf(fp, '%f\t',lamda_curvelet);
    fprintf(fp, '\n');
        fprintf(fp,'%3s ', 'RLtime=');
    fprintf(fp, '%f\t',RLtime);
    fprintf(fp, '\n');
end

%% Comment on parameter choice
fprintf(fp,'%3s ', '****Reminder of parameter choice****');
fprintf(fp, '\n');
paraflag=0;
if flag1==1
    if lambda1>0.3
            fprintf(fp,'%3s ', '****Reminder of parameter choice****');
            fprintf(fp,'The input 3D-framelet thresholding value is too large, please select a lower one.');
            fprintf(fp, '\n');
            fprintf(fp,'If lower parameters cannot effectively reduce noise, please carefully examine the image, as there may be artifacts formed by noise that has not been well removed.');
            fprintf(fp, '\n');
            paraflag=1;           
    end
end
if flag2==1&&flag1~=1
    if sparsity_curvelet>0.75
        if fun_all(1)<fun_all(2)
             fprintf(fp,'Noise is reduced in the image, but the blurring is not alleviated.');
            fprintf(fp, '\n');
            fprintf(fp,'The input image SNR seems high, a lower framelet&curvelet thresholding value could be chosen to further improve resolution.');
            fprintf(fp, '\n');
            paraflag=1;
        end
    end
    if sparsity_curvelet<0.7
        if fun_all(2)<0.1*fun_all(1)
            fprintf(fp,'The input image SNR seems not sufficiently high, please take care of the possible artifacts caused by residual noise in the image.');
            fprintf(fp, '\n');
             fprintf(fp,'It is recommended to try a higher framelet&curvelet thresholding value.');
            fprintf(fp, '\n');
            paraflag=1;
        end
    end
    if max(lamda_framelet,lamda_curvelet)<0.1*761.5*exp(-17.24*sparsity_curvelet)-0.0001584
            fprintf(fp,'The input image SNR seems not sufficiently high, please take care of the possible artifacts caused by residual noise in the image.');
            fprintf(fp, '\n');
             fprintf(fp,'It is recommended to try a higher framelet&curvelet thresholding value.');
            fprintf(fp, '\n');
            paraflag=1;
    end
end
if paraflag==0
        fprintf(fp,'No warning about the parameter choice.');
        fprintf(fp, '\n');
        fprintf(fp,'If you want to further improve resolution, try lowering the framelet&curvelet thresholding.');
        fprintf(fp, '\n');
        fprintf(fp,'On the contrary, if you want to further reduce noise, try increasing the framelet&curvelet thresholding.');
        fprintf(fp, '\n'); 
else
            msgbox('Possible inappropriate parameters or extreme situations detected, please check the comments in the .tif file accompanying the deconvolution results','Warning')
end
fclose(fp);




% --- Executes on button press in Activation2D.
function Activation2D_Callback(hObject, eventdata, handles)
% hObject    handle to Activation2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Activation2D



function lambda_framelet_Callback(hObject, eventdata, handles)
% hObject    handle to lambda_framelet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lambda_framelet as text
%        str2double(get(hObject,'String')) returns contents of lambda_framelet as a double


% --- Executes during object creation, after setting all properties.
function lambda_framelet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lambda_framelet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lambda_curvelet_Callback(hObject, eventdata, handles)
% hObject    handle to lambda_curvelet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lambda_curvelet as text
%        str2double(get(hObject,'String')) returns contents of lambda_curvelet as a double


% --- Executes during object creation, after setting all properties.
function lambda_curvelet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lambda_curvelet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RLtimes_Callback(hObject, eventdata, handles)
% hObject    handle to RLtimes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RLtimes as text
%        str2double(get(hObject,'String')) returns contents of RLtimes as a double


% --- Executes during object creation, after setting all properties.
function RLtimes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RLtimes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Activation3D.
function Activation3D_Callback(hObject, eventdata, handles)
% hObject    handle to Activation3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Activation3D



function lambda1_Callback(hObject, eventdata, handles)
% hObject    handle to lambda1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lambda1 as text
%        str2double(get(hObject,'String')) returns contents of lambda1 as a double


% --- Executes during object creation, after setting all properties.
function lambda1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lambda1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxframe_Callback(hObject, eventdata, handles)
% hObject    handle to maxframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxframe as text
%        str2double(get(hObject,'String')) returns contents of maxframe as a double


% --- Executes during object creation, after setting all properties.
function maxframe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Wavelength_Callback(hObject, eventdata, handles)
% hObject    handle to Wavelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Wavelength as text
%        str2double(get(hObject,'String')) returns contents of Wavelength as a double


% --- Executes during object creation, after setting all properties.
function Wavelength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Wavelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Pixelsize_Callback(hObject, eventdata, handles)
% hObject    handle to Pixelsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pixelsize as text
%        str2double(get(hObject,'String')) returns contents of Pixelsize as a double


% --- Executes during object creation, after setting all properties.
function Pixelsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pixelsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NA_Callback(hObject, eventdata, handles)
% hObject    handle to NA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NA as text
%        str2double(get(hObject,'String')) returns contents of NA as a double


% --- Executes during object creation, after setting all properties.
function NA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Scale_Callback(hObject, eventdata, handles)
% hObject    handle to Scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Scale as text
%        str2double(get(hObject,'String')) returns contents of Scale as a double


% --- Executes during object creation, after setting all properties.
function Scale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global filename;
global pathname;
global bit;
global f_stack_old;

%% Check if MRA-related operation is activated
flag1=get(handles.Activation3D,'value');
flag2=get(handles.Activation2D,'value');
if flag2~=1&&flag1~=1
     msgbox('Auto MRA para is used for MRA parameter initilization, please  select at least one MRA option','Warning')
end
if exist('f_stack_old')==1
f=f_stack_old(:,:,1);
[~,~,d3]=size(f_stack_old);
f=single(f);
f=f/max(max(f));
[~,sparsity_curvelet] = Sparsitycal_fun(f);
sparsity_curvelet
if flag1==1&&d3==1
     msgbox('Single frame input,only 2D MRA parameters are automatically decided','Warning')
end
if flag1==1&&d3>1
    a=20;b =11.26;c=4.352e-05;predicted= a*exp(-b*sparsity_curvelet)+c;
    [predicted] = Rounding(predicted);
        if mean(mean(f))<0.1
            predicted=predicted*0.5;
        end
    handles.lambda1.String=num2str(predicted);
   if flag2==1
               handles.lambda_framelet.String=num2str('0.005');
               handles.lambda_curvelet.String=num2str('0.005');
   end
elseif flag2==1
    a = 761.5;b =17.24; c = -0.0001584;predicted= a*exp(-b*sparsity_curvelet)+c;
    if predicted<1e-6
        predicted=1e-6;
    end
    [predicted] = Rounding(predicted);
    handles.lambda_framelet.String=num2str(predicted);
    handles.lambda_curvelet.String=num2str(predicted);
end
else 
    msgbox('No image selected')
end
