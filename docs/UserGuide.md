# User guide

## Introduction

FRETzel performs cell segmentation and FRET analysis of cells over a time series in a MATLAB guided user interface (GUI). Useful for identifying cells from brightfield images in messy and crowded data sets, such as those produced from the growth of adipocytes (fat cells), but can also be used on other cells, such as yeast. 

See diplayGUI.jpg and completeGUI.jpg for the user interface at opening and once an experiment using adipocytes has been analysed. 


## Using the GUI

Ensure all the code is downloaded (see Getting started) and added to MATLAB path using 'set path'. The process of the GUI is as follows:

#### 1. Loading files
- Run the FRETzel.m and the GUI will open (you should see displayGUI.jpg). Instructions are on screen. The user clicks 'Get file' to select the first file in the time series or the pre frame. The brightfield channel and FRET channels are displayed on screen. 
-At this point the user should check the parameters are correct (see Parameters below). 

#### 2. Cell selection, segmentation and analysis
- The user clicks 'Select cells'. A cursor appears and the user selects desired cells in the brightfield image with a left click, backspace undoes the previous selected point and right click on the final cell ends the cell selection. 
- Upon right clicking, circles (ellipses) appear around each cell and can be changed in shape and size by the user to best fit the outline of each cell. Zoom in, zoom out and grab (in the top left corner) can be used at this point to move around the brightield image if necessary (the image needs to be zoomed back out once complete). For some non-uniform cells, it can be difficult to find a perfect fit, but the cell segmentation process can generally cope with this. Once happy, the user presses enter.
- The GUI finds cell outlines and plots them on the brightfield image, the cells are numbered and have their radius in microns. 
- The user clicks 'Run analysis'. Graphs of FRET ratio, donor intensity, FRET intensity and acceptor intensity for each cell with time (in minutes) and a cell size histogram are plotted. See completeGUI.jpg.
- FRET ratio is calculated as FRET intensity/donor intensity for each cell at each time frame. 

#### 3. Saving files
- The user clicks 'Save as' to save 5 files in a selected folder. The GUI will overwrite previously saved files of the same name if they are in the selected folder. cellImage.png, FRETgraphs.png, cellSizeHistogram.png, outputFRETdata.mat and outputFRETdata.xlsx. MATLAB figures open and close themselves during saving.
- cellImage.png contains the brightfield image, labelled with the cell outlines and radii, the donor image, the FRET image and acceptor image of the first file with titles. 
- FRETgraphs.png contains graphs of FRET ratio, donor intensity, FRET intensity and acceptor intensity with time (in minutes).
- cellSizeHistogram.png contains the cell size histogram.
- outputFRETdata.mat saves a MATLAB workspace containing variables 'maskCellFinal' (cell outlines for each cell), 'FRETInt' (the FRET channel intensity for each cell at each time frame), 'donorInt' (the donor channel intensity for each cell at each time frame), 'acceptorInt' (the acceptor channel intensity for each cell at each time frame), 'timesMins' (the time of each frame in minutes) and 'cellRadii' (the radius of each cell outlined). These can then be used for further data analysis.
- outputFRETdata.xlsx saves an Excel worksheet with cell radii and FRET ratio, FRET intensity, donor intensity, and acceptor intensity for each cell at each time frame. 

## Calculating FRET data

The software uses the time stamps from each of the .lsm files in the user selected folder to order the data. The final cell outlines are used for each of the fluorescence channels in each file to plot intensity of each fluorescence (FRET, donor and acceptor) and the FRET ratio (FRET fluorescence/donor fluorescence). 

## Parameters

See top right corner of displayGUI.jpg for parameters needed to run the GUI. Descriptions of each are below.

#### File type
Preset to 'lsm', can be any of the Bio-formats file formats. Must be input without a dot infront ('lsm' not '.lsm').

#### Approximate diameter of cells
This must be input in microns. It is preset to 20 microns. It is used to calculate the size of the ellipses created in the cell segmentation process, based on the pixel size of the image.

#### Active contour number of iterations
This is a number between 0 and 100, preset to 50 that determines the number of iterations MATLAB function `activecontour` carries out. `activecontour` segments the cells once the user has confirmed the ellipses to be the desired size and shape. For messier data sets (i.e those using adipocytes or that have noise) a lower number of iterations (around 20) is more useful for identifying cell outlines. For data sets with low noise and more uniform cells (i.e. yeast cells) a higher number of iterations (70-100) is more useful. This number can be changed to find the optimum for each data set.

#### Number of channels
The number of channels is the total number of channels in each file, which can be more than the 4 specified in 'File requirements' above, but needs to be specified for the analysis of the cells across all files. Preset to 4.

#### Input time interval or read time interval from files
There is the option for the user to set a uniform time interval between channels (in seconds), or to tick the box 'Read time interval from files' to get the time data from the metaData of each file. One of these options must be input to plot the donor, FRET and acceptor graphs with time. 

If the user inputs a time interval, then it is important that the files in the folder have uniform names that have consecutive numbering in order of when they were aquired. For example - control01.lsm, control02.lsm, ..., control10.lsm, control11.lsm (but not control1.lsm, control2.lsm, ..., control10.lsm, control11.lsm).

#### First frame is a 'pre' frame
If this box is ticked, the GUI will set the first file in the time series to time = -5 minutes, assuming it is a 'pre' frame (i.e. before the addition of of any solution used in the time series to measure the change in FRET). It sets the second file to time = 0 minutes. If it is not ticked, the GUI will set the first file to time = 0 minutes. In both cases each consecutive file is set using the time interval specified or found in the metaData. 

## File requirements
Software can use any Bio-Formats file formats (i.e. .tiff, .tif, .lsm). Channels within the files must be in a certain order: donor fluorescence, FRET fluorescence, acceptor fluorescence then a brightfield/scattered light image. There can be more than these 4 channels but the FRET and brightfield channels must be the first 4 and in the correct order. For each time series the files must be in a separate folder.  

The GUI can process four different types and combinations of time series:

1. Separate files for each time frame of the time series (for example, control01.lsm, control02.lsm, control03.lsm... in a single folder). Where there is a pre frame. 
2. Separate files for each time frame of the time series (for example, control01.lsm, control02.lsm, control03.lsm... in a single folder). Where there is not a pre frame. 
3. A time series combined in one file (for example, a folder that contains control01.lsm, that has the time frames and respective channels in one file).
4. A separate pre frame file and a time series combined in one file (for example, a folder that contains control01.lsm as the pre frame and control02.lsm, that has the remaining time frames and respective channels in one file).

A pre frame is a frame taken before the addition of of any solution used in the time series to measure the change in FRET, it must have the same field of view as the following files in the time series. The GUI can differentiate between these 4 combinations using the user set parameters (see above).

The GUI can be slower with files 50000KB and above, but still works. 

## Information about downloads required

### FRETzel.m
Code used to run the GUI, made in MATLAB's GUIDE. Loads Bio-formats files using ImEx1.m, segments cells using `imellipse`, `createMask` and `activecontour` and runs the analysis of cells, calculating the FRET ratio and plotting graphs.

GUIDE will be removed in a future release of MATLAB. The software will still be usable but not editable. At some point the GUI will be migrated or made an executable file. 

### FRETzel.fig
The GUI window.

### bfmatlab 
MATLAB Bio-Formats toolbox from Open Microscopy Environment. 

### ImEx1.m
Uses Bio-Formats toolbox to open and access data from microscopy files.
