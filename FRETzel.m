function varargout = FRETzel(varargin)
% FRETZEL MATLAB code for FRETzel.fig
%      FRETZEL, by itself, creates a new FRETZEL or raises the existing
%      singleton*.
%
%      H = FRETZEL returns the handle to a new FRETZEL or the handle to
%      the existing singleton*.
%
%      FRETZEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FRETZEL.M with the given input arguments.
%
%      FRETZEL('Property','Value',...) creates a new FRETZEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FRETzel_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FRETzel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FRETzel

% Last Modified by GUIDE v2.5 07-May-2021 15:09:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FRETzel_OpeningFcn, ...
                   'gui_OutputFcn',  @FRETzel_OutputFcn, ...
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

% --- Executes just before FRETzel is made visible.
function FRETzel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FRETzel (see VARARGIN)

% Choose default command line output for FRETzel
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = FRETzel_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% Creates following text, edit boxes and axes.
function axes1_CreateFcn(hObject, eventdata, handles)   

function text3_CreateFcn(hObject, eventdata, handles)   % Instructions

function text10_CreateFcn(hObject, eventdata, handles)  % Parameters

function text8_CreateFcn(hObject, eventdata, handles)   % Approx cell diameter

function text9_CreateFcn(hObject, eventdata, handles)   % Active contour number of iterations

function text7_CreateFcn(hObject, eventdata, handles)   % Acceptor channel

function text6_CreateFcn(hObject, eventdata, handles)   % FRET channel

function text5_CreateFcn(hObject, eventdata, handles)   % Donor channel

function text12_CreateFcn(hObject, eventdata, handles)  % File:

function edit1_Callback(hObject, eventdata, handles)    % approx cell size in microns

function edit1_CreateFcn(hObject, eventdata, handles)   % approx cell size in microns

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit2_Callback(hObject, eventdata, handles)    % active contour number of iterations

function edit2_CreateFcn(hObject, eventdata, handles)    % active contour number of iterations

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit3_Callback(hObject, eventdata, handles)     % number of channels

function edit3_CreateFcn(hObject, eventdata, handles)   % number of channels

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function checkbox1_Callback(hObject, eventdata, handles)  % first frame is 'pre' frame.
% get(handles.checkbox1,'Value') returns toggle state of checkbox1 (ie. 0
% if unticked, 1 if ticked.

function checkbox3_Callback(hObject, eventdata, handles)        % Read in time from files


function edit4_Callback(hObject, eventdata, handles)            % File type (i.e. .tiff or .lsm)

function edit4_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit5_Callback(hObject, eventdata, handles)            % Time interval in seconds

function edit5_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1 - Get File. User has to select 
% 1st file from selected folder. Each of the 4 channels in the file are placed 
% into one of the axes created above with brightfield on the main axes.
function pushbutton1_Callback(hObject, eventdata, handles)
handles = guidata(hObject);                     % opens handles structure

cla(handles.axes1, 'reset')                     % resets
handles.axes1.Visible = 'off'; 
cla(handles.axes2)
cla(handles.axes3)
cla(handles.axes4)
cla(handles.axes6)
cla(handles.axes7)
cla(handles.axes8)
cla(handles.axes9)
cla(handles.axes11)
handles.text3.Visible = 'on'; 
handles.text6.Visible = 'off';
handles.text5.Visible = 'off';
handles.text7.Visible = 'off';
handles.text12.Visible = 'off';

set(handles.text3, 'String', 'Check file type. Click "Get file" to begin.')

[FileName, folder]= uigetfile('*');         % user selects first file in folder of choice

if FileName == 0                            % Prevents error if user presses cancel

else
    set(handles.text3, 'String', 'Loading file...')
if ispc==1
    [imageData,metaData]=imEx1(strcat(folder,'\',FileName));
else
        [imageData,metaData]=imEx1(strcat(folder,'/',FileName));
end
    handles.imageData = imageData;                  
    handles.folder = folder;
    handles.metaData = metaData;
    imageDataN=imageData;                           

    brightfield=cat(3,mat2gray(imageDataN(:,:,4)),mat2gray(imageDataN(:,:,4)),mat2gray(imageDataN(:,:,4)));     % saves individual channels
    handles.brightfield = brightfield;
    donorImage=cat(3,zeros(size(imageDataN(:,:,4))),mat2gray(imageDataN(:,:,1)),mat2gray(imageDataN(:,:,1)));
    handles.donorImage = donorImage;
    FRETImage=cat(3,mat2gray(imageDataN(:,:,2)),mat2gray(imageDataN(:,:,2)),zeros(size(imageDataN(:,:,4))));
    handles.FRETImage = FRETImage;
    acceptorImage=cat(3,mat2gray(imageDataN(:,:,3)),mat2gray(imageDataN(:,:,3)),zeros(size(imageDataN(:,:,4))));
    handles.acceptorImage = acceptorImage;

    axes(handles.axes1);                               % Plots 4 channels on different axes
    check4 = get(handles.checkbox4, 'Value'); 
if check4==0
    imshow(brightfield, 'Parent', handles.axes1)
    handles.im=brightfield;
else
    imshow(brightfield+donorImage+FRETImage, 'Parent', handles.axes1)
    handles.im=brightfield+donorImage+FRETImage;
end
    
    handles.text5.Visible = 'on';
    axes(handles.axes2);
    imshow(donorImage, 'Parent', handles.axes2)
    handles.text6.Visible = 'on';
    axes(handles.axes3)
    imshow(FRETImage, 'Parent', handles.axes3)
    handles.text7.Visible = 'on';
    axes(handles.axes4)
    imshow(acceptorImage, 'Parent', handles.axes4)
    drawnow

    myString = sprintf('File: %s%s', folder, FileName);
    handles.text12.String = myString;
    handles.text12.Visible = 'on';
    set(handles.text3, 'String', 'Set approximate cell size above. Click "Select cells"')
    pixelSize=double(metaData.getPixelsPhysicalSizeX(0).value());           % pixel size from meta data
sizeUnit=char(metaData.getPixelsPhysicalSizeX(0).unit().getSymbol());   % units from data - should be in microns
handles.sizeUnit = sizeUnit;
handles.pixelSize = pixelSize;
handles.cellRadius = []; % initialise here to add later
handles.maskCellFinal = [];
end

guidata(hObject, handles);

% --- Executes on button press in pushbutton9 - select cells.
function pushbutton9_Callback(hObject, eventdata, handles)

handles = guidata(hObject); 
imageData = handles.imageData;
metaData = handles.metaData;
pixelSize=handles.pixelSize;
sizeUnit=handles.sizeUnit;
set(handles.text3, 'String', 'Left click each cell to select. Click backspace to undo previous cell selection. Right click final cell to finish.')

[x, y] = getpts(handles.axes1);         % cell selection

cellApproxDiam = str2double(get(handles.edit1,'string'));               % takes input from user of cell diameter in microns
% pixelSize=double(metaData.getPixelsPhysicalSizeX(0).value());           % pixel size from meta data
% sizeUnit=char(metaData.getPixelsPhysicalSizeX(0).unit().getSymbol());   % units from data - should be in microns 
sizeThresh=100; 
cellDiam = cellApproxDiam/pixelSize;

for  p=1:length(x)                    % calculates ellipses around cells and allows them to be changed in size/shape. Ends when user presses enter
    Spot_coords=[x(p), y(p)];
    hE(p)=imellipse(gca,[Spot_coords(1)-cellDiam/2,Spot_coords(2)-cellDiam/2,cellDiam,cellDiam]);       %ellipse size based on user input cell diameter
end
set(handles.text3, 'String', 'Adjust each circle size and shape to outline each cell. Use zoom and grab to move image. Press enter when finished.')
drawnow
%pause
uiwait(msgbox('Press ok when finished adjusting outlines'));
set(handles.text3, 'String', 'Drawing cell outlines...')
drawnow
niterations = str2double(get(handles.edit2,'string'));                      

for p=1:length(x)             % finds and plots approximated cell outline 
   
    maskCell=createMask(hE(p));
    bw = activecontour(imageData(:,:,4),maskCell(1:size(imageData(:,:,4),1),1:size(imageData(:,:,4),2)),niterations);       %number of iterations of active contour depends on user input
    bw2 = imfill(bw(:,:,1),'holes');    % fills in holes in image by flooding
    se = strel('disk',2);               % disk structuring element
    bw3 = imopen(bw2, se );             % morphological opening to remove objects smaller than structuring element
    bw4 = bwareaopen(bw3, sizeThresh);
    cc = bwconncomp(bw4,4);
    L_cells = labelmatrix(cc);
    sizeMask=[];
    for l=1:max(L_cells(:))
        sizeMask(l)=sum(L_cells(:)==l);
    end
    [~,I]=max(sizeMask);
    bw5=L_cells==I;
    [row,col]=find(bwperim(bw5));
    hold on
    scatter(col,row,5,'.')
    maskCellFinal(:,:,p)=bw5(1:size(imageData,1),1:size(imageData,2));
    rp=regionprops(maskCellFinal(:,:,p),'Area','Centroid');
    cellRadius(p)=(rp.Area/pi)^0.5;
    text(rp.Centroid(1),rp.Centroid(2),strcat('c',num2str(p),':',num2str(cellRadius(p)*pixelSize,3),sizeUnit),'Color','White')
end
delete(hE);                 % Deletes ellipses

set(handles.text3, 'String', 'Check parameters above. Then click "Run analysis"')
oldCellRadius=handles.cellRadius;
updatedCellRadius=cat(2,oldCellRadius,cellRadius);
handles.cellRadius = updatedCellRadius;
oldMask=handles.maskCellFinal;
updatedMask=cat(3, oldMask,maskCellFinal);
%size(updatedMask)
handles.maskCellFinal = updatedMask;
% handles.sizeUnit = sizeUnit;
% handles.pixelSize = pixelSize;
guidata(hObject, handles);


% --- Executes on button press in pushbutton10 - run analysis.
function pushbutton10_Callback(hObject, eventdata, handles)
cla(handles.axes6)
cla(handles.axes7)
cla(handles.axes8)
cla(handles.axes9)
cla(handles.axes11)
handles = guidata(hObject); 
maskCellFinal = handles.maskCellFinal;
folder = handles.folder;
imageData1 = handles.imageData;
cellRadius = handles.cellRadius;
pixelSize = handles.pixelSize;
metaData = handles.metaData;

set(handles.text3, 'String', 'Analysing data...')
drawnow

check3 = get(handles.checkbox3, 'Value');           % finds value of 'read time interval from file' tickbox
check1 = get(handles.checkbox1, 'Value');           % finds value of pre frame tickbox
fileType = get(handles.edit4,'string');             % Gets file type
fileType1 = strcat('*.',fileType);                  
timeInt = str2double(get(handles.edit5,'string'));  % Time interval in seconds
imageFiles=dir(fullfile(folder, fileType1));        % all image files from previously selected folder in struct variable imageFiles
numChannels = str2double(get(handles.edit3,'string'));  % Gets number of channels from user

if check3 == 0                                                              % uses the user input time interval in seconds
    
    if check1 == 1                                                          % there is a pre frame
if ispc==1
        [imageData2,~]=imEx1(strcat(folder,'\',imageFiles(2).name));        % reads second file in folder
else
        [imageData2,~]=imEx1(strcat(folder,'/',imageFiles(2).name));        % reads second file in folder
end
        if size(imageData2, 3) > numChannels                                % pre frame and timeseries    

            imageDataN = cat(3, imageData1, imageData2);                    % concatenates pre frame and timeseries 
            numFiles = size(imageData2, 3);                                 
            times = -300;                                                   % sets pre-frame time at -5 min (in seconds) 
            timesMins = [times, (0:timeInt:(numFiles/numChannels-1)*timeInt)]/60;   % sets times for each frame in minutes 
            for j = 1:numFiles/numChannels+1                                % for the number of frames that there are in total
                for p=1:size(maskCellFinal,3)                               % for the number of cell outlines
                    donorImage=imageDataN(:,:,1+(j-1)*numChannels);          % donor channel in time frame
                    FRETimage=imageDataN(:,:,2+(j-1)*numChannels);          % FRET channel in time frame
                    acceptorImage=imageDataN(:,:,3+(j-1)*numChannels);        % acceptor channel in time frame 
                    donorInt(p,j)=sum(donorImage(maskCellFinal(:,:,p)==1));     % donor intensity of cell
                    FRETInt(p,j)=sum(FRETimage(maskCellFinal(:,:,p)==1));      % FRET intensity of cell
                    acceptorInt(p,j)=sum(acceptorImage(maskCellFinal(:,:,p)==1));  % acceptor intensity of cell
                end
            end 
            handles.timesMins = timesMins'; 

        elseif size(imageData2, 3) == numChannels                           % pre frame and no timeseries 

            for j=1:length(imageFiles)
                if ispc==1
                [imageData,~]=imEx1(strcat(folder,'\',imageFiles(j).name)); 
                else
                [imageData,~]=imEx1(strcat(folder,'/',imageFiles(j).name)); 
                end
                
                
                for p=1:size(maskCellFinal,3)
                    donorImage=imageData(:,:,1);
                    FRETimage=imageData(:,:,2);
                    acceptorImage=imageData(:,:,3);
                    donorInt(p,j)=sum(donorImage(maskCellFinal(:,:,p)==1));
                    FRETInt(p,j)=sum(FRETimage(maskCellFinal(:,:,p)==1));
                    acceptorInt(p,j)=sum(acceptorImage(maskCellFinal(:,:,p)==1));
                end
            end
            times = -300;                  
            timesMins = [times, (1:timeInt:(length(imageFiles)-1)*timeInt)]/60;
            timesMins(2) = 0;
            handles.timesMins = timesMins';   

        end

    elseif check1 == 0                                                      % there isn't a pre frame

        numFiles = size(imageData1, 3);
        if numFiles > numChannels                                           % no pre frame and timeseries - only works for one time series

            for j = 1:numFiles/numChannels    
                for p=1:size(maskCellFinal,3)
                    donorImage=imageData1(:,:,1+(j-1)*numChannels);
                    FRETimage=imageData1(:,:,2+(j-1)*numChannels);
                    acceptorImage=imageData1(:,:,3+(j-1)*numChannels);
                    donorInt(p,j)=sum(donorImage(maskCellFinal(:,:,p)==1));
                    FRETInt(p,j)=sum(FRETimage(maskCellFinal(:,:,p)==1));
                    acceptorInt(p,j)=sum(acceptorImage(maskCellFinal(:,:,p)==1));
                end
            end
            timesMins = (0:timeInt:(numFiles/numChannels-1)*timeInt)/60;
            handles.timesMins = timesMins'; 

        elseif numFiles == numChannels                                      % no pre frame and no timeseries

            for j=1:length(imageFiles)
                if ispc==1
                [imageData,~]=imEx1(strcat(folder,'\',imageFiles(j).name));
                else
                [imageData,~]=imEx1(strcat(folder,'/',imageFiles(j).name));
                end
                for p=1:size(maskCellFinal,3)
                    donorImage=imageData(:,:,1);
                    FRETimage=imageData(:,:,2);
                    acceptorImage=imageData(:,:,3);
                    donorInt(p,j)=sum(donorImage(maskCellFinal(:,:,p)==1));
                    FRETInt(p,j)=sum(FRETimage(maskCellFinal(:,:,p)==1));
                    acceptorInt(p,j)=sum(acceptorImage(maskCellFinal(:,:,p)==1));
                end
            end
            timesMins = (0:timeInt:(length(imageFiles)-1)*timeInt)/60;
            handles.timesMins = timesMins';  
        end 
 
    end

else                                                                        % use time from metaData instead of user input time interval
    
    if check1 == 1                                                          % there is a preframe

        if ispc==1
        [imageData2,metaData2]=imEx1(strcat(folder,'\',imageFiles(2).name));
        else
                    [imageData2,metaData2]=imEx1(strcat(folder,'/',imageFiles(2).name));
        end
        if size(imageData2, 3) > numChannels                                % pre frame and timeseries    
            times = -300;                                                   
            imageDataN = cat(3, imageData1, imageData2);
            numChannels = str2double(get(handles.edit3,'string')); 
            numFiles = size(imageData2, 3);

            for i = 0:numChannels:numFiles-numChannels
                times=[times, double(metaData2.getPlaneDeltaT(0,i).value())];   % time data from timeseries metaData
            end
            timesMins = times/60;

            for j = 1:numFiles/numChannels+1   
                for p=1:size(maskCellFinal,3)
                    donorImage=imageDataN(:,:,1+(j-1)*numChannels);
                    FRETimage=imageDataN(:,:,2+(j-1)*numChannels);
                    acceptorImage=imageDataN(:,:,3+(j-1)*numChannels);
                    donorInt(p,j)=sum(donorImage(maskCellFinal(:,:,p)==1));
                    FRETInt(p,j)=sum(FRETimage(maskCellFinal(:,:,p)==1));
                    acceptorInt(p,j)=sum(acceptorImage(maskCellFinal(:,:,p)==1));
                end
            end 
            handles.timesMins = timesMins'; 

        elseif size(imageData2, 3) == numChannels                            % pre frame and no timeseries 

            for j=1:length(imageFiles)
                if ispc==1
                [imageData,metaData]=imEx1(strcat(folder,'\',imageFiles(j).name));
                else
                [imageData,metaData]=imEx1(strcat(folder,'/',imageFiles(j).name)); 
                end
                timeString=metaData.getImageAcquisitionDate(0).string();    % gets time data from each individual file infolder 
                t(j) = datetime(timeString,'InputFormat','yyyy-MM-dd''T''HH:mm:ss');
                for p=1:size(maskCellFinal,3)
                    donorImage=imageData(:,:,1);
                    FRETimage=imageData(:,:,2);
                    acceptorImage=imageData(:,:,3);
                    donorInt(p,j)=sum(donorImage(maskCellFinal(:,:,p)==1));
                    FRETInt(p,j)=sum(FRETimage(maskCellFinal(:,:,p)==1));
                    acceptorInt(p,j)=sum(acceptorImage(maskCellFinal(:,:,p)==1));
                end
            end
            [~,Ind]=sort(t);                                                % sorts t into ascending order
            handles.t = t;                  
            t2=t(Ind);
            handles.t2 = t2;                
            imageFiles=imageFiles(Ind);                                     % orders files into ascending order
            timesMins=zeros(length(imageFiles),1);
            timesMins(1)=-5;                                                % sets pre frame time
            timesMins(2)=0;                                                 % sets first time to 0
            timesMins(3:length(imageFiles))=minutes(t2(3:end)-t2(2));  

            donorInt=donorInt(:,Ind);                                       % orders variables
            FRETInt=FRETInt(:,Ind);
            acceptorInt=acceptorInt(:,Ind);
            handles.timesMins = timesMins; 
        end

    elseif check1 == 0                                                      % there is no pre frame

        numFiles = size(imageData1, 3);
        if numFiles > numChannels                                           % no pre frame and timeseries - only works for one time series
            times = [];
            for i = 0:numChannels:numFiles-numChannels
                times=[times, double(metaData.getPlaneDeltaT(0,i).value())];
            end
            timesMins = times/60;
            for j = 1:numFiles/numChannels    
                for p=1:size(maskCellFinal,3)
                    donorImage=imageData1(:,:,1+(j-1)*numChannels);
                    FRETimage=imageData1(:,:,2+(j-1)*numChannels);
                    acceptorImage=imageData1(:,:,3+(j-1)*numChannels);
                    donorInt(p,j)=sum(donorImage(maskCellFinal(:,:,p)==1));
                    FRETInt(p,j)=sum(FRETimage(maskCellFinal(:,:,p)==1));
                    acceptorInt(p,j)=sum(acceptorImage(maskCellFinal(:,:,p)==1));
                end
            end
            handles.timesMins = timesMins'; 

        elseif numFiles == numChannels                                      % no pre frame and no timeseries

            for j=1:length(imageFiles)
                if ispc==1
                [imageData,metaData]=imEx1(strcat(folder,'\',imageFiles(j).name));
                else
              [imageData,metaData]=imEx1(strcat(folder,'/',imageFiles(j).name));
                end
                timeString=metaData.getImageAcquisitionDate(0).string();
                t(j) = datetime(timeString,'InputFormat','yyyy-MM-dd''T''HH:mm:ss');
                for p=1:size(maskCellFinal,3)
                    donorImage=imageData(:,:,1);
                    FRETimage=imageData(:,:,2);
                    acceptorImage=imageData(:,:,3);
                    
                    
                    donorInt(p,j)=sum(donorImage(maskCellFinal(:,:,p)==1));
                    FRETInt(p,j)=sum(FRETimage(maskCellFinal(:,:,p)==1));
                    acceptorInt(p,j)=sum(acceptorImage(maskCellFinal(:,:,p)==1));
                end
            end
            [~,Ind]=sort(t);                
            handles.t = t;                  
            t2=t(Ind);
            handles.t2 = t2;                
            imageFiles=imageFiles(Ind);     
            timesMins=zeros(length(imageFiles),1);
            timesMins(1)=0;                    
            timesMins(2:length(imageFiles))=minutes(t2(2:end)-t2(1));  

            donorInt=donorInt(:,Ind);       
            FRETInt=FRETInt(:,Ind);
            acceptorInt=acceptorInt(:,Ind);
            handles.timesMins = timesMins;  
            
        end         
    end       
end

axes(handles.axes6)             % plots FRET ratio
for p=1:size(maskCellFinal,3)
    plot(timesMins,FRETInt(p,:)./donorInt(p,:))
    xlabel('Time (minutes)')
    ylabel('FRET ratio')
    hold on
end 
errorbar(timesMins, mean(FRETInt./donorInt,1),std(FRETInt./donorInt,0,1)/size(maskCellFinal,3)^0.5,'k','LineWidth',1)

axes(handles.axes7)             % plots donor intensity 
for p=1:size(maskCellFinal,3)
    plot(timesMins,donorInt(p,:))
    xlabel('Time (minutes)')
    ylabel('Donor intensity')
    hold on
end 
errorbar(timesMins, mean(donorInt,1),std(donorInt,0,1)/size(maskCellFinal,3)^0.5,'k','LineWidth',1)

axes(handles.axes8)             % plots FRET intensity
for p=1:size(maskCellFinal,3)
    plot(timesMins,FRETInt(p,:))
    xlabel('Time (minutes)')
    ylabel('FRET intensity')
    hold on
end
errorbar(timesMins, mean(FRETInt,1),std(FRETInt,0,1)/size(maskCellFinal,3)^0.5,'k','LineWidth',1)


axes(handles.axes9)             % plots acceptor intensity
for p=1:size(maskCellFinal,3)
    plot(timesMins,acceptorInt(p,:))
    xlabel('Time (minutes)')
    ylabel('Acceptor intensity')
    hold on
end
errorbar(timesMins, mean(acceptorInt,1),std(acceptorInt,0,1)/size(maskCellFinal,3)^0.5,'k','LineWidth',1)


axes(handles.axes11)            % plots cell size histogram
cellRadii = cellRadius*pixelSize;
histogram(cellRadii)
xlabel('Cell Radii (\mum)')
ylabel('Count')

set(handles.text3, 'String', 'Click "Save as" to select folder to save graphs and data')

handles.donorInt = donorInt;     
handles.FRETInt = FRETInt;
handles.acceptorInt = acceptorInt;

guidata(hObject, handles);


% --- Executes on button press in pushbutton7 - Save as. User chooses
% directory  to save data in and CellImage.png, FRETgraphs.png,
% cellSizeHistogram.png, outputFRETdata.mat and outputFRETdata.xlsx are saved.
function pushbutton7_Callback(hObject, eventdata, handles)

handles = guidata(hObject);                 
brightfield = handles.brightfield;
donorImage = handles.donorImage;
FRETImage = handles.FRETImage;
acceptorImage = handles.acceptorImage;
maskCellFinal = handles.maskCellFinal;
acceptorInt = handles.acceptorInt;
FRETInt = handles.FRETInt;
donorInt = handles.donorInt;
timesMins = handles.timesMins;
sizeUnit = handles.sizeUnit;
pixelSize = handles.pixelSize;

folder = uigetdir;                                                          % User chooses directory to save data in

set(handles.text3, 'String', 'Saving cell images...')
CellImage = figure('Renderer', 'painters', 'Position', [10 10 1500 1000]);  % Makes a new figure to save CellImage.png
subplot(2,2,1)
imshow(brightfield)
title('Brightfield')
hold on
for p=1:size(maskCellFinal,3)                                               % Adds cell outlines and annotations to brightfield plot
    [row,col]=find(bwperim(maskCellFinal(:,:,p)));
    scatter(col,row,5,'.')
    rp=regionprops(maskCellFinal(:,:,p),'Area','Centroid');
    cellRadius(p)=(rp.Area/pi)^0.5; 
    text(rp.Centroid(1),rp.Centroid(2),strcat('c',num2str(p),':',num2str(cellRadius(p)*pixelSize,3),sizeUnit),'Color','White')
    hold on
end
subplot(2,2,2)
imshow(donorImage)
title('Donor')
subplot(2,2,3)
imshow(FRETImage)
title('FRET')
subplot(2,2,4)
imshow(acceptorImage)
title('Acceptor')
saveas(CellImage, fullfile(folder,'CellImage.fig'))
close(CellImage)   

set(handles.text3, 'String', 'Saving graphs...')
FRETgraphs = figure;                                                        % Makes a new figure to save FRETgraph.png
for p=1:size(maskCellFinal,3)
    subplot(2,2,1)
    plot(timesMins,FRETInt(p,:)./donorInt(p,:))
    xlabel('Time (minutes)')
    ylabel('FRET ratio')
    hold on

    subplot(2,2,2)
    plot(timesMins,donorInt(p,:))
    xlabel('Time (minutes)')
    ylabel('Donor intensity')
    hold on
    
    subplot(2,2,3)
    plot(timesMins,FRETInt(p,:))
    xlabel('Time (minutes)')
    ylabel('FRET intensity')
    hold on
    
    subplot(2,2,4)
    plot(timesMins,acceptorInt(p,:))
    xlabel('Time (minutes)')
    ylabel('Acceptor intensity')
    hold on
end
saveas(FRETgraphs, fullfile(folder,'FRETgraphs.fig'))
      
close(FRETgraphs)

cellRadii = cellRadius*pixelSize;
cellSizeHistogram = figure;
histogram(cellRadii)
xlabel('Cell Radii (\mum)')
ylabel('Count')
saveas(cellSizeHistogram, fullfile(folder,'cellSizeHistogram.png'))         % Makes a new figure to save cellSizeHistogram.png
close(cellSizeHistogram)

set(handles.text3, 'String', 'Saving data in MATLAB workspace...')
save(fullfile(folder, 'outputFRETdata.mat'),'maskCellFinal','FRETInt','donorInt', 'acceptorInt', 'timesMins', 'cellRadii') % Saves outputFRETdata.mat 

set(handles.text3, 'String', 'Saving data in excel file...')
filename = fullfile(folder, 'outputFRETdata.xlsx');                         % Saves outputFRETdata.xlsx

xlswrite(filename,1:size(maskCellFinal,3),'cellRadii','1');
xlswrite(filename,cellRadius*pixelSize,'cellRadii','2');

finalRatio=(FRETInt./donorInt);
xlswrite(filename,timesMins,'FRETratio','B');
xlswrite(filename,finalRatio','FRETratio','C');

xlswrite(filename,timesMins,'donor','B');
xlswrite(filename,donorInt','donor','C');

xlswrite(filename,timesMins,'FRETint','B');
xlswrite(filename,FRETInt','FRETint','C');

xlswrite(filename,timesMins,'acceptor','B');
xlswrite(filename,acceptorInt','acceptor','C');

set(handles.text3, 'String', 'Done. Press "Reset" or "Get file" to start again')


% --- Executes on button press in pushbutton8 - Reset. Axes are cleared and
% user can begin process again.
function pushbutton8_Callback(hObject, eventdata, handles) 
cla(handles.axes1, 'reset')
handles.axes1.Visible = 'off'; 
cla(handles.axes2)
cla(handles.axes3)
cla(handles.axes4)
cla(handles.axes6)
cla(handles.axes7)
cla(handles.axes8)
cla(handles.axes9)
cla(handles.axes11)

set(handles.text3, 'String', 'Click "Get file" to begin')
handles.text6.Visible = 'off';
handles.text5.Visible = 'off';
handles.text7.Visible = 'off';
handles.text12.Visible = 'off';


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('1')
brightfield = handles.brightfield;
donorImage = handles.donorImage;
FRETImage = handles.FRETImage;
acceptorImage = handles.acceptorImage;
allFimage=mat2gray(donorImage+FRETImage+acceptorImage);

range=get(hObject,'Max')-get(hObject,'Min');
(get(hObject,'Value')-get(hObject,'Min'))/range;
% mask=thresholdSegmentLabel(mat2gray(handles.im),50,1,(get(hObject,'Value')-get(hObject,'Min'))/range);
threshVal=graythresh(allFimage);

cellApproxDiam = str2double(get(handles.edit1,'string'));               % takes input from user of cell diameter in microns
% pixelSize=double(metaData.getPixelsPhysicalSizeX(0).value());           % pixel size from meta data
% sizeUnit=char(metaData.getPixelsPhysicalSizeX(0).unit().getSymbol());   % units from data - should be in microns 

cellAreaLim = 0.1*(cellApproxDiam/(2*handles.pixelSize))^2*pi;

mask=thresholdSegmentLabel(allFimage,round(cellAreaLim),1,((get(hObject,'Value')-get(hObject,'Min'))/range)*threshVal);

handles.mask=mask;
%     [col,row]=find(bwperim(mask>0));
%     axes(handles.axes1); 
%     hold off
%     imshow(handles.im,[])
%     hold on
%     scatter(row,col,'.')

    axes(handles.axes1); 
    hold off
    imshow(handles.im,[])
    hold on
   
for l=1:max(mask(:))
    [col,row]=find(bwperim(mask==l));
    scatter(row,col,1,'y')
end
guidata(hObject, handles);

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles) %select thresholded cells
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text3, 'String', 'Left click each cell to select. Click backspace to undo previous cell selection. Right click final cell to finish.')

imageData = handles.imageData;
metaData = handles.metaData;
pixelSize=handles.pixelSize;
sizeUnit=handles.sizeUnit;


[x, y] = getpts(handles.axes1); 
set(handles.text3, 'String', 'Separating individual cells')

seedMask=zeros(size(handles.im,1),size(handles.im,2));
for f=1:length(x)
    seedMask(round(y(f)),round(x(f)))=1;
end
%seedMask=imdilate(seedMask,strel('disk',5));
mask=handles.mask;
mask=mask>0;
%NewMask=WatershedCellLabels(handles.mask,mask,seedMask);
waterImage=~mask-seedMask;
waterImage(~mask)=Inf;
NewMask = watershed(waterImage);
NewMask(~mask) = 0;

size(NewMask)

    axes(handles.axes1); 
% for l=1:max(NewMask(:))
%     [col,row]=find(bwperim(NewMask==l));
%     scatter(row,col,1)
% end
p=1;
cellRadius=[];
    maskCellFinal=zeros(size(handles.im,1),size(handles.im,2),length(x));
    for lbl = 1:max(NewMask(:))  
        sum(sum((NewMask ==lbl).*seedMask))>0
%p=lbl;
if sum(sum((NewMask ==lbl).*seedMask))>0
            maskCellFinal(:,:,p) = NewMask == lbl; %# find pixels belonging to current label. First object is the borders I think, so take the objects after these as indivudal cell masks  
    rp=regionprops(maskCellFinal(:,:,p),'Area','Centroid');
    cellRadius(p)=(rp.Area/pi)^0.5;
      [col,row]=find(bwperim(maskCellFinal(:,:,p)));
    scatter(row,col,1)
    text(rp.Centroid(1),rp.Centroid(2),strcat('c',num2str(p),':',num2str(cellRadius(p)*pixelSize,3),sizeUnit),'Color','White')
    p=p+1;
end
    end

%cellApproxDiam = str2double(get(handles.edit1,'string'));               % takes input from user of cell diameter in microns
% pixelSize=double(metaData.getPixelsPhysicalSizeX(0).value());           % pixel size from meta data
% sizeUnit=char(metaData.getPixelsPhysicalSizeX(0).unit().getSymbol());   % units from data - should be in microns 

set(handles.text3, 'String', 'Check parameters above. Then click "Run analysis"')
oldCellRadius=handles.cellRadius;
updatedCellRadius=cat(2,oldCellRadius,cellRadius);
handles.cellRadius = updatedCellRadius;
oldMask=handles.maskCellFinal;
updatedMask=cat(3, oldMask,maskCellFinal);
%size(updatedMask)
handles.maskCellFinal = updatedMask;
% handles.sizeUnit = sizeUnit;
% handles.pixelSize = pixelSize;
guidata(hObject, handles);


