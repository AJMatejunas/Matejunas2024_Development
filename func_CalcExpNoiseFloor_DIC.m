function [NF] = func_CalcExpNoiseFloor_DIC(StaticPath,StaticFile,BitD,grid,specimen, ...
    time,material,shearExtrap,shearSmooth,bulkExtrap,bulkSmooth,diffOpts,...
    globalOpts,gridMethodOpts)
%This script is written to calculate the Noise floor for the grey-levels 
    %and kinematic fields: Diplacement, strain, acceleration, and stress 
    %gauge stress 
 
%Author: Andrew Matejunas

%Version History/changelog
    %2023-07-06: Original version

%Function Input Arguments;
    %StaticPath
    %BitD- Number of bits encoding the image grey-levels (generally 8 or
        %16)
   %Grid- Structure containing grid method parameters
   %specimen- structure containing specimen geometry
   %time- standard structure containing time information
   %material- structure containing known material properties
   %ShearExtrap- extrapolation options used in the full identification
        %procedure for the shear portion of the cost function
   %ShearSmooth- Smoothing options used in the full identifiation procedure
        %for the full shear modulus identification procedure
   %BulkExtrap- Extrapolation options for the bulk modulus identification
   %BulkSmooth- Smoothing options for the bulk modulus identification
   %diffOpts- differentiation settings to determine acceleration and strain
        %rate
   %globalOpts- global processing options
   %gridMethodOpts
    

%function outputs:
    %NF- structure containg a record of experimental noise floors with
        %fields:
            %GLA- Grey level noise floor (absolute)
            %GL- Grey level noise (percentage of dynamic range)
            %Shear.disp.x- displacement in the x
            %Shear.disp.y- displacement in the y
            %disp.Bulk.x- displacement in the x
            %disp.Bulk.y- displacement in the y
            %Shear.strain.x- axial strain in the x direction
            %Shear.strain.xAv- average strain in the x
            %Shear.strain.y- transverse strain in the y direction
            %Shear.strain.yAv- average strain in the y
            %Shear.strain.s- in-plane shear strains
            %Shear.strain.sAv- average in-plane shear stress
            %Bulk.strain.x- axial strain in the x direction
            %Bulk.strain.xAv- average strain in the x
            %Bulk.strain.y- transverse strain in the y direction
            %Bulk.strain.yAv- average strain in the y
            %Bulk.strain.s- in-plane shear strains
            %Bulk.strain.sAv- average in-plane shear stress
            %Shear.accel.x- acceleration in the x direction
            %Shear.accel.y- acceleration in the y direction
            %Bulk.accel.x- acceleration in the x direction
            %Bulk.accel.y- acceleration in the y direction
            %Shear.SG.x- stress gage stress in the x direction using shear
                %identification parameters
            %Bulk.SG.x- stress gage stress in the x direction using bulk
                %identification parameters
            %Shear.SG.s- stress gage stress in in-plane shear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Determine number of grey levels in the images
TotGL=2^BitD;

%% Load images
imageNoise.addNoise=false;

%% Calculate Grey Level Noise floor
fprintf('Calculating GrayLevel Noise Level \n')

%Load Reference image
refIm=imread([StaticPath,StaticFile]);

%Determine repeating part of the image string
[~,~,imageExt] = fileparts(StaticFile);
 usLoc = find(StaticFile == '_'); %find underscores
% %Repeating part is everything before the final underscore
startStr = StaticFile(1:usLoc(end)); 

% Get the names of all image files with the same ext in the directory
fstart = startStr;
StaticFileStruct = dir([StaticPath,fstart,'*',imageExt]);
StaticFileCell = {StaticFileStruct.name}';

% Sort the data files in numerical order
sortedFiles=cell(size(StaticFileCell));
sortedFiles{1} = StaticFile;
for i = 2:length(StaticFileCell)
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

    for j = 2:length(StaticFileCell)       
        if strcmp(checkStr1,StaticFileCell{j}) || strcmp(checkStr2,StaticFileCell{j})
            sortedFiles{i} = StaticFileCell{j};
            break
        end
    end
end
StaticFileCell = sortedFiles;


%Read Second image
ref2=imread([StaticPath,StaticFileCell{2}]);

%Calculate Grey Level Noise floor;
imDiff=double(refIm-ref2);
NF.GLA=std(imDiff,0,'all'); %absolute graylevel noise floor
NF.GL=NF.GLA/TotGL*100; %percentage of overall dynamic range


%% Calculate displacement noise floor
fprintf('Calculating Displacement Noise Floor \n')
% fprintf('GRID METHOD PROCESSING\n')
%fprintf('--------------------------------------------------------------\n')
%--------------------------------------------------------------------------
% GRID IMAGE PROCESSING

% fprintf('Processing images using the grid method toolbox.\n')
% Process the image squence with the grid method toolbox
[grid,pos,disp] = func_gridMethodImageProcessing(StaticPath,...
    StaticFile,...
    grid,gridMethodOpts,imageNoise);

%--------------------------------------------------------------------------
% Update Geometry and Number of Frames Based on Displacement Matrix Size
%fprintf('Updating geometry variables based on the size of displacement field matrices.\n')
[specimen,~] = func_updateSpecGeom(specimen,grid,disp);

% Currently the rotations are unused so remove them to save RAM
%disp = rmfield(disp,'rot');


%--------------------------------------------------------------------------
% Create the time vector based on the number of frames in the disp struct
% time.numFrames = size(disp.x,3);
%time.vec = 0:time.step:(size(disp.x,3)-1)*time.step;

%Create arrays of x and y vectors
% X_vec=pos.x;
% Y_vec=pos.y;
pos.lengthX = pos.x(end)+pos.xStep/2;
pos.lengthY = pos.y(end)+pos.yStep/2;
pos.xGridF = padarray(pos.xGrid,[0,0,time.numFrames-1], ...
    'replicate','post');
pos.yGridF = padarray(pos.yGrid,[0,0,time.numFrames-1], ...
    'replicate','post');
pos.x0F = squeeze(padarray(pos.x,[0,0,time.numFrames-1], ...
    'replicate','post'));

fprintf('Obtaining and setting the free edge location.\n')
[~,~,disp] =...
    func_getFreeEdge(globalOpts.hardCodeFreeEdge,...
    StaticPath,StaticFile,specimen,disp);

% Extrapolate the data to account for the missing pitch on the edges
% Crop and Extrapolate displacement fields

%Noise Floor using Shear Processing parameters
fprintf('Croping and Extrapolating Displacement Fields \n')
[Shear.disp.x,Shear.disp.y,Shear.disp.rX,Shear.disp.rY]=...
    func_cropAndExtrapDispFields_v4(pos,disp.x,disp.y, ...
    shearExtrap.disp,false);

NF.Shear.disp.x=std(Shear.disp.x(:,:,2),0,"all");
NF.Shear.disp.y=std(Shear.disp.y(:,:,2),0,"all");

%Displacement Noise Floor using Axial Processing Parameters
[Bulk.disp.x,Bulk.disp.y,Bulk.disp.rX,Bulk.disp.rY]=...
    func_cropAndExtrapDispFields_v4(pos,disp.x,disp.y, ...
    bulkExtrap.disp,false);

NF.Bulk.disp.x=std(Bulk.disp.x(:,:,2),0,"all");
NF.Bulk.disp.y=std(Bulk.disp.y(:,:,2),0,"all");



%% Deterimine Strain Noise Floors
fprintf('Determining Strain Noise Floor \n')
%For Shear
[Shear.strain,Shear.disp]=func_smoothCalcStrain_v4(pos,time,Shear.disp,...
    shearSmooth.strain,shearExtrap.strain,false);

NF.Shear.strain.x=std(Shear.strain.x(:,:,2),0,'all');
NF.Shear.strain.xAvg=std(Shear.strain.xAvg(:,2),0,'all');
NF.Shear.strain.y=std(Shear.strain.y(:,:,2),0,"all");
NF.Shear.strain.yAvg=std(Shear.strain.yAvg(:,2),0,"all");
NF.Shear.strain.s=std(Shear.strain.s(:,:,2),0,"all");
NF.Shear.strain.sAvg=std(Shear.strain.sAvg(:,2),0,"all");


%For Bulk
[Bulk.strain,Bulk.disp]=func_smoothCalcStrain_v4(pos,time,Bulk.disp,...
    bulkSmooth.strain,bulkExtrap.strain,false);

NF.Bulk.strain.x=std(Bulk.strain.x(:,:,2),0,'all');
NF.Bulk.strain.xAvg=std(Bulk.strain.xAvg(:,2),0,'all');
NF.Bulk.strain.y=std(Bulk.strain.y(:,:,2),0,"all");
NF.Bulk.strain.yAvg=std(Bulk.strain.yAvg(:,2),0,"all");
NF.Bulk.strain.s=std(Bulk.strain.s(:,:,2),0,"all");
NF.Bulk.strain.sAvg=std(Bulk.strain.sAvg(:,2),0,"all");

%% Determine Acceleration 
fprintf('Determining Acceleration Noise Floor \n')
[Shear.accel,~,Shear.disp] = func_smoothCalcAccel_v4(pos,time,Shear.disp, ...
    shearSmooth.accel,shearExtrap.accel,...
    diffOpts,false);

NF.Shear.accel.x=std(Shear.accel.x(:,:,2),0,'all');
NF.Shear.accel.y=std(Shear.accel.y(:,:,2),0,'all');

[Bulk.accel,~,Bulk.disp] = func_smoothCalcAccel_v4(pos,time,Bulk.disp, ...
    bulkSmooth.accel,bulkExtrap.accel,...
    diffOpts,false);

NF.Bulk.accel.x=std(Bulk.accel.x(:,:,2),0,'all');
NF.Bulk.accel.y=std(Bulk.accel.y(:,:,2),0,'all');

%% Determine Stress-Gauge Noise Floor 
fprintf('Calculating the stress gauge Noise Floor \n')
Shear.SG=func_Full_SG(Shear.accel,pos.x,time,material.rho);
Bulk.SG=func_Full_SG(Bulk.accel,pos.x,time,material.rho);

NF.Shear.SG.x=std(Shear.SG.x(:,2),0,"all");
NF.Shear.SG.x=std(Shear.SG.s(:,2),0,"all");

NF.Bulk.SG.x=std(Bulk.SG.x(:,2),0,"all");
NF.Bulk.SG.x=std(Bulk.SG.s(:,2),0,"all");



end