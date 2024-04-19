function [grid,pos,disp] = func_gridMethodImageProcessing_AJM(refImagePath,refImageFile,...
    grid,gridMethodOpts,imageNoise)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton

% Adapted by: Andrew Matejunas

% Date: 9/3/2017
%ChangeLog:
    %2022-12-05: Locked in use of an existing specimen location file to
        %specify image region of interest

% Processes a sequence of grid images to obtain displacement fields
% This code uses the grid processing tool box that can be found at:
% www.thegridMethodOpts.net, developed by Grediac

% Don't add noise to images by default
if nargin < 5
    imageNoise.addNoise = false; % if there's not enough inputs, don't add noise
end

%--------------------------------------------------------------------------
% 1) Load the Reference Image and Find All Images in Folder
% Get the repeated part of the string assuming underscore termination
[~,~,imageExt] = fileparts(refImageFile);   % determine image file extension, .tif
usLoc = find(refImageFile == '_'); % stores index location of underscore in image file name
if ~isempty(usLoc)  % if it finds an underscore, starter string is from beginning to the underscore
    startStr = refImageFile(1:usLoc(end));
else
    startStr = 'DefGrid_';
end

% Ask for the filename starter and new filename
% AJM Commenting this part out to allow for uninterrupted running of my
% noise sweep code

% resp = inputdlg({'Enter the repeating part of the image file string:'},...
%              'Input repeating file string', 1, {startStr} );  
% needs to know repeating portion of image file names      
% fstart = resp{1}; % stores start of file name
fstart=startStr;

% Get the names of all image files with the same ext in the directory
refImageFileStruct = dir([refImagePath,fstart,'*',imageExt]); % '*' means any set of characters in between
refImageFileCell = {refImageFileStruct.name}'; % returns names of files, stores in cell array

% Sort the data files in numerical order
sortedFiles{1} = refImageFile;  % places the first image file as first entry in sorted image cell array
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
% 2) Select specimen ROI to mask images
%%fprintf('Getting specimen region of interest (ROI).\n')

specLocFile = 'specimenLocation.mat';
%AJM- commenting this part out to automatically load existing ROI file

% loadSpecLocFile = false;
% % Check if the specimen location file already exists, if it does prompt to load it
% if exist([refImagePath,specLocFile],'file') == 2
%     choice = questdlg('Specimen location (ROI) file found:', ...
%     'Process Grid Images', 'Load ROI file','Reselect ROI','Load ROI file');
%     switch choice
%         case 'Load ROI file'
%             loadSpecLocFile = true;
%         case 'Reselect ROI'
%             loadSpecLocFile = false;
%     end
% end
    %%fprintf('Loading existing ROI file.\n')
    load([refImagePath,specLocFile])


% Assign the position of the specimen to the pos struct for return to main
pos.specimenLoc.bottomLeft = specimenLoc.bottomLeft;
pos.specimenLoc.topRight = specimenLoc.topRight;

%--------------------------------------------------------------------------
% 3) Input Grid Parameters
%AJM- Maing grid parameters automatic
%%fprintf('Input grid parameters.\n')
% gridData = inputdlg({'Number of pixels per period:','Grid pitch: (m)'}, ...
%              'Grid analysis parameters', 1, {num2str(grid.pxPerPeriod), num2str(grid.pitch)} );
% grid.pxPerPeriod = str2double(gridData{1}); 
% grid.pitch = str2double(gridData{2});
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
%%fprintf('Spatial unwrapping.\n')
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
    range.x = (specimenLoc.bottomLeft(1):specimenLoc.topRight(1));
    range.y = (specimenLoc.bottomLeft(2):specimenLoc.topRight(2));
    phi.x_nuw = phi.x;
    phi.y_nuw = phi.y;
    
    %%fprintf('Temporal unwrappping.\n')
    threshold = pi;
    phase = func_temporalUnwrap(phi,threshold,'field',range);

    for i = 1:size(phi.x,3)
        phi.x(:,:,i) = phi.x(:,:,i) + 2*pi*phase.x(i); 
        phi.y(:,:,i) = phi.y(:,:,i) + 2*pi*phase.y(i);
    end
end

if gridMethodOpts.debug
    figure;
    hold on
    plot(squeeze(mean(mean(phi.x(range.y,range.x,:)))),'-+b')
    plot(squeeze(mean(mean(phi.x_nuw(range.y,range.x,:)))),'-xr')
    xlabel('Frame')
    ylabel('Mean Phase Over ROI')
    legend('Temp Unwrapped','No Unwrap')
    hold off    
end

% Calculate Displacements and Strains
%%fprintf('Calculating displacement and strain components.\n')
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
        squeeze(phi.x(dispRangeY,dispRangeX,i)),squeeze(phi.y(dispRangeY,dispRangeX,i)),gridMethodOpts.dispCalcMethod,100);  
end

%--------------------------------------------------------------------------
% 6) Convert Displacement and Crop to ROI 
% Convert the displacement to mm
%fprintf('Cropping to ROI and updating position/size vectors.\n')
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
    % Update grid size parameters
    grid.length = size(disp.x,2)*grid.pitchX/grid.pxPerPeriod;
    grid.height = size(disp.y,1)*grid.pitchY/grid.pxPerPeriod; 
    
    % Create the position mesh grid
    grid = rmfield(grid,'mPerPx');
    grid.mPerPxX = grid.pitchX/grid.pxPerPeriod;
    grid.mPerPxY = grid.pitchY/grid.pxPerPeriod;
    pos.x = grid.mPerPxX/2:grid.mPerPxX:(grid.mPerPxX*size(disp.x,2));
    pos.y = grid.mPerPxY/2:grid.mPerPxY:(grid.mPerPxY*size(disp.x,1));
    %pos.y = flip(pos.y - grid.height/2);
    [pos.xGrid,pos.yGrid] = meshgrid(pos.x,pos.y);
    pos.xStep = grid.mPerPxX;
    pos.yStep = grid.mPerPxY;

else
    % Update grid size parameters
    grid.length = size(disp.x,2)*grid.pitch/grid.pxPerPeriod;
    grid.height = size(disp.y,1)*grid.pitch/grid.pxPerPeriod;
    
    % Create the position mesh grid
    pos.x = grid.mPerPx/2:grid.mPerPx:(grid.mPerPx*size(disp.x,2));
    pos.y = grid.mPerPx/2:grid.mPerPx:(grid.mPerPx*size(disp.x,1));
    %pos.y = flip(pos.y - grid.height/2);
    [pos.xGrid,pos.yGrid] = meshgrid(pos.x,pos.y);
    pos.xStep = grid.mPerPx;
    pos.yStep = grid.mPerPx;    
end

%--------------------------------------------------------------------------
% Ask user about grid rotation and correct if needed
% Grid Rotation is not utilized commenting out rotation correction
% %fprintf('Correcting for any grid rotation.\n')
% [grid,disp] = func_rotateDisp(grid,disp,~(gridMethodOpts.hardCodeRotAngle));

%--------------------------------------------------------------------------

%FOLLOWING CODE REMOVED FOR COMPUTATIONAL EFFICIENCY
% Tranpose fields if the images are in portrait
% if isfield(grid,'transposeFields')
%     if grid.transposeFields
%         temp = disp.x;
%         disp.x = permute(disp.y,[2,1,3]);
%         disp.y = permute(temp,[2,1,3]);
%         temp = pos.xGrid;
%         pos.xGrid = permute(pos.yGrid,[2,1,3]);
%         pos.yGrid = permute(temp,[2,1,3]);
%         temp = pos.x;
%         pos.x = pos.y;
%         pos.y = temp;
%     end
% end

end

