function [grid,pos,disp] = func_gridMethodProcessingImageDef(refImagePath,refImageFile,...
    grid,gridMethodOpts,specLoc,imageNoise)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 17/1/2017
% Date Edited: 4/5/2019
%
% Processes a sequence of grid images to obtain displacement fields
% This code uses the grid processing tool box that can be found at:
% www.thegridMethodOpts.net, developed by Grediac

%--------------------------------------------------------------------------
% 1) Load the Reference Image and Find All Images in Folder
% Get the repeated part of the string assuming underscore termination
[~,~,imageExt] = fileparts(refImageFile);
usLoc = find(refImageFile == '_');
if ~isempty(usLoc)
    startStr = refImageFile(1:usLoc(end));
else
    startStr = 'DefGrid_';
end       
fstart = startStr;

% Get the names of all image files with the same ext in the directory
refImageFileStruct = dir([refImagePath,fstart,'*',imageExt]);
refImageFileCell = {refImageFileStruct.name}';

% Sort the data files in numerical order
sortedFiles{1} = refImageFile;
for i = 2:length(refImageFileCell)
    checkStr1 = [fstart,num2str(i),imageExt];   
    
    nn = length(num2str(i));
    if nn==1
        num = ['00' num2str(i)];
    elseif nn==2
        num = ['0' num2str(i)];
    else
        num = num2str(i);
    end
    checkStr2 = [fstart,num,imageExt];
    
    for j = 2:length(refImageFileCell)       
        if strcmp(checkStr1,refImageFileCell{j}) || strcmp(checkStr2,refImageFileCell{j})
            sortedFiles{i} = refImageFileCell{j};
            break
        end
    end
end
refImageFileCell = sortedFiles;

% Load the reference image 
refImage = imread([refImagePath,refImageFile]);
refImage = double(refImage);

if imageNoise.addNoise
    refImage = func_addNoiseToImagesStruct(refImage,imageNoise);
end
    
%--------------------------------------------------------------------------
% 2) Set range for image masking
range.x = (specLoc.bottomLeft(1):specLoc.topRight(1));
range.y = (specLoc.bottomLeft(2):specLoc.topRight(2));
grid.mPerPx = grid.pitch/grid.pxPerPeriod;

%--------------------------------------------------------------------------
% 4) Build the Analysis Window
%Build Window
analysisWindow = build_window(gridMethodOpts.windowFlag,...
    gridMethodOpts.windowWidth*grid.pxPerPeriod); % See documentation
% 0 = Gaussian window, default
% 1 = Bi-triangular window, localised phenom, low noise

%--------------------------------------------------------------------------
% 5) Process Images
% Pre-alloc vars for speed
[sy,sx] = size(refImage);
st = length(refImageFileCell );
disp.x = zeros([sy,sx,st]); 
disp.y = zeros([sy,sx,st]);
disp.rot = zeros([sy,sx,st]);
phi.x = zeros([sy,sx,st]);
phi.y = zeros([sy,sx,st]);

% Spatial Unwrapping
for i = 1:length(refImageFileCell)
    % Load the current image file
    currImage = imread([refImagePath,refImageFileCell{i}]);
    currImage = double(currImage);
    
    % Add noise to the image if required
    if imageNoise.addNoise
        currImage = func_addNoiseToImagesStruct(currImage,imageNoise);
    end

    % Calculation of the phase and spatial phase unwrapping
    [phi.x(:,:,i),phi.y(:,:,i),phaseMod.x(:,:,i),phaseMod.y(:,:,i)] = ...
        LSA(currImage, analysisWindow, grid.pxPerPeriod);
    phi.x(:,:,i) = unwrap2D(single(phi.x(:,:,i)));
    phi.y(:,:,i) = unwrap2D(single(phi.y(:,:,i)));
end

% Temporal Unwrapping
% CAUTION: this only works for a sequence of images with small increments in
% the phase, uses average phase over the FOV to calc rigid body motion
if gridMethodOpts.temporalUnwrap
    range.x = (specLoc.bottomLeft(1):specLoc.topRight(1));
    range.y = (specLoc.bottomLeft(2):specLoc.topRight(2));
    phi.x_nuw = phi.x;
    phi.y_nuw = phi.y;
    
    threshold = pi;
    phase = func_temporalUnwrap(phi,threshold,'field',range);

    for i = 1:size(phi.x,3)
        phi.x(:,:,i) = phi.x(:,:,i) + 2*pi*phase.x(i); 
        phi.y(:,:,i) = phi.y(:,:,i) + 2*pi*phase.y(i);
    end
end

% Calculate Displacements and Strains
disp.x = zeros(size(phi.x));
disp.y = zeros(size(phi.x));
disp.rot = zeros(size(phi.x));

% Calculate a window over which the specimen exists + an extra pitch to
% avoid nans - currently uses the range variable
% X Range
if (min(range.x)-grid.pxPerPeriod) < 1
    startInd = 1;
else
    startInd = min(range.x)-grid.pxPerPeriod;
end
if (max(range.x)+grid.pxPerPeriod) > size(phi.x,2)
    endInd = size(phi.x,2);
else
    endInd = max(range.x)+grid.pxPerPeriod;
end
dispRangeX = startInd:endInd;

% Y Range
if (min(range.y)-grid.pxPerPeriod) < 1
    startInd = 1;
else
    startInd = min(range.y)-grid.pxPerPeriod;
end
if (max(range.y)+grid.pxPerPeriod) > size(phi.y,1)
    endInd = size(phi.y,1);
else
    endInd = max(range.y)+grid.pxPerPeriod;
end
dispRangeY = startInd:endInd;

for i = 2:length(refImageFileCell)  
    [disp.x(dispRangeY,dispRangeX,i), disp.y(dispRangeY,dispRangeX,i),~,~,~, disp.rot(dispRangeY,dispRangeX,i)]...
        = calculate_U_EPS(grid.pxPerPeriod,squeeze(phi.x(dispRangeY,dispRangeX,1)),squeeze(phi.y(dispRangeY,dispRangeX,1)),...
        squeeze(phi.x(dispRangeY,dispRangeX,i)),squeeze(phi.y(dispRangeY,dispRangeX,i)),gridMethodOpts.dispCalcMethod,100,true);  
end

%--------------------------------------------------------------------------
% 6) Convert Displacement and Crop to ROI 
% Convert the displacement to mm
if grid.asymmPitch
    disp.x = disp.x*(grid.pitchX/grid.pxPerPeriod); 
    disp.y = disp.y*(grid.pitchY/grid.pxPerPeriod);
else
    disp.x = disp.x*(grid.pitch/grid.pxPerPeriod); 
    disp.y = disp.y*(grid.pitch/grid.pxPerPeriod);
end

% Crop the image to the ROI
disp.x = disp.x(range.y,range.x,:);
disp.y = disp.y(range.y,range.x,:);

if grid.asymmPitch
    % Create the position mesh grid
    grid = rmfield(grid,'mPerPx');
    grid.mPerPxX = grid.pitchX/grid.pxPerPeriod;
    grid.mPerPxY = grid.pitchY/grid.pxPerPeriod;
    pos.x = grid.mPerPxX/2:grid.mPerPxX:(grid.mPerPxX*size(disp.x,2));
    pos.y = grid.mPerPxY/2:grid.mPerPxY:(grid.mPerPxY*size(disp.x,1));
    [pos.xGrid,pos.yGrid] = meshgrid(pos.x,pos.y);
    pos.xStep = grid.mPerPxX;
    pos.yStep = grid.mPerPxY;

    % Update grid size parameters
    grid.length = size(disp.x,2)*grid.pitchX/grid.pxPerPeriod;
    grid.height = size(disp.y,1)*grid.pitchY/grid.pxPerPeriod; 
else
    % Create the position mesh grid
    pos.x = grid.mPerPx/2:grid.mPerPx:(grid.mPerPx*size(disp.x,2));
    pos.y = grid.mPerPx/2:grid.mPerPx:(grid.mPerPx*size(disp.x,1));
    [pos.xGrid,pos.yGrid] = meshgrid(pos.x,pos.y);
    pos.xStep = grid.mPerPx;
    pos.yStep = grid.mPerPx;

    % Update grid size parameters
    grid.length = size(disp.x,2)*grid.pitch/grid.pxPerPeriod;
    grid.height = size(disp.y,1)*grid.pitch/grid.pxPerPeriod;
end

%--------------------------------------------------------------------------
% Ask user about grid rotation and correct if needed
% Suppress user input, hard code rotation angle
[grid,disp] = func_rotateDisp(grid,disp,false);

end

