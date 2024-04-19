function [OutputStruct] =...
func_InterpFEtoGM(FEDataFile,SaveDir,Desig,subset)
%This function is written to export FE displacements interpolated onto grid
    %Method coordinates.

%Original Author: Lloyd Fletcher
%Adapted by: Andrew Matejunas

%Date completed: 
%Changelog/version history:
    %2022/10/02- added subset input. For this function it would typically
        %be 1
        %made camera size=output displacement fields from the grid method
            %processing algorithm
        %Changed to a lagrangian reference frame

%Function input arguments:
    %FEDataFile- string containing the file name, with path, for the File
        %containing finite element displacement data in the original finite
        %element coordinate system
    %SaveDir- Directory to save function outputs
    %Desig- Designation for the simulations
    %subset-% Pixel subsampling for the FE and grid data, default of 3 is fine for most
        %applications. Large values will slow computation drastically.

    


%Function output arguments:
    %OutputStruct- Structure containing the original FE deformations and
        %the interpolated displacements
    %


%% Load FE data file
fprintf('Loading Image Deformation data \n')
load(FEDataFile);

%% INITIALISE: Data Structures for Camera Parameters
fprintf('Intialising data structure for the camera and the grid.\n')

%//////////////////////////////////////////////////////////////////////////
% USER SPECIFIED VALUES:
% Parameters for the simulated camera
cameraData.px_x = 389;              % Number of pixels on the sensor in the x direction
cameraData.px_y = 246;              % Number of pixels on the sensor in the y direction
cameraData.bit = 10;                % Number of bits encoding the image HPVX uses 10
cameraData.contrastOffset = 0.5;    % Fraction of dynamic range for centre of grey level distrib    
cameraData.contrastAmp = 0.325;     % Fraction of dynamic range, grey level amplitude of grid pattern  
% NOTE: the contrast amplitude above is typical of an appropriately blurred
% image taken using the HPVX

% Parameters for the grid
grid.pxPerPeriod = 5;   % in pixels per period
grid.pitch = 0.9e-3;    % grd pitch in meters

% Specimen location within the FOV
FEData.specLocX =0;     % X location of the bottom left hand corner in the FOV in m
FEData.specLocY = 0;   % Y location of the bottom left hand corner in the FOV in m
%//////////////////////////////////////////////////////////////////////////

% Pixel subsampling for the FE and grid data, default of 3 is fine for most
% applications. Large values will slow computation drastically.


% Create remaining grid parameters based on the FEData structure
grid.mPerPx = grid.pitch/grid.pxPerPeriod;
grid.length = FEData.Lx;
grid.width = FEData.Ly;
grid.specLocX = FEData.specLocX;
grid.specLocY = FEData.specLocY;
grid.specLocXPx = ceil(FEData.specLocX/grid.mPerPx);
grid.specLocYPx = ceil(FEData.specLocY/grid.mPerPx);

% Create remaining camera parameters
cameraData.numImages = FEData.numFrames;
cameraData.dynRange = 2^cameraData.bit;
numLinesX = cameraData.px_x/grid.pxPerPeriod;
numLinesY = cameraData.px_y/grid.pxPerPeriod;
cameraData.FOV_x = numLinesX*grid.pitch;
cameraData.FOV_y = numLinesY*grid.pitch;
cameraData.px_res = cameraData.FOV_y/cameraData.px_y;

%% Convert FE data to Camera Sensor Size
fprintf('Interpolating FE displacement data onto pixel array:\n')
% Create mesh of nodal locations
X = (0:FEData.elemXSize:(FEData.Lx+FEData.elemXSize/2))+FEData.specLocX;
Y = (0:FEData.elemYSize:(FEData.Ly+FEData.elemYSize/2))+FEData.specLocY;
[FEData.Xn,FEData.Yn] = meshgrid(X,Y);
clear X Y

% Create mesh of pixel edge values (when upsampling edges are not lost)
X = cameraData.px_res/(2*subset):cameraData.px_res/subset:...
    cameraData.FOV_x-cameraData.px_res/(2*subset);
Y = cameraData.px_res/(2*subset):cameraData.px_res/subset:...
    cameraData.FOV_y-cameraData.px_res/(2*subset);
[cameraData.subPxLocXm,cameraData.subPxLocYm] = meshgrid(X,Y);
clear X Y

% Interpolate FE data to be the same size as the camera pixel array
% including any subsampling required
for i = 1:FEData.numFrames
    fprintf('Interpolating FE displacement frame: %i\n',i)
    x = FEData.Xn;
    y = FEData.Yn;
    
    ux = squeeze(FEData.dispX(:,:,i));
    uy = squeeze(FEData.dispY(:,:,i));
    
    scatInterpUx = scatteredInterpolant(x(:),y(:),ux(:),'natural','none');
    scatInterpUy = scatteredInterpolant(x(:),y(:),uy(:),'natural','none');
    FEData.subPxDispX(:,:,i) = scatInterpUx(cameraData.subPxLocXm,cameraData.subPxLocYm);
    FEData.subPxDispY(:,:,i) = scatInterpUy(cameraData.subPxLocXm,cameraData.subPxLocYm);
end
FEData.disp.x=FEData.dispX;
FEData.disp.y=FEData.dispY;

%% Extract displacements to expected data structuee
disp.x=FEData.subPxDispX;
disp.y=FEData.subPxDispY;



%% Save
SaveFile=strcat(SaveDir,'/',Desig,'_FEdispInterp2Grid.mat');
save(SaveFile)

%% Set result data to the propper output structure
OutputStruct.File=SaveFile;
OutputStruct.disp=disp;
OutputStruct.grid=grid;
OutputStruct.cameraData=cameraData;
OutputStruct.FEData=FEData;
OutputStruct.x=x;
OutputStruct.y=y;


end