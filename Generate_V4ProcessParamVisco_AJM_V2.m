%This script is written to generate processing parameters for viscoelastic
%IBII experiments. It has been optimized to feed into V4 of the main IBII
%codes particularly func_smoothCalcStrain_V4 and func_SmoothCalcAccel_V4

%Author:  Andrew Matejunas
%Date created: 2023/02/28

%VersionHistory/ChangeLog
    %Initial Version: Generates processing parameters used to calculate
        %displacement, accelerations, strains, and strain rates from grid 
        %images
    %2023/03/05- Added inputs for the data conditioning before passing into
        %the cost function
                %Added SolveOpts data structure to determine how the cost
                    %function is minimized
                 %Added minOpts structure to govern how fmincon is run
   %V2-2023/08/03- Adapted to work with 
              %IBII_ViscoStandardSolid_ExperimentalKG_IdentV2
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize
clear variables; close all; clc

%% Choose save location
SaveDir=uigetdir('','Choose location to save processing parameters');

%% Input test designation
procDeg=char(cell2mat(inputdlg('input test designation')));

%% Deterimine Source of the displacement fields
    %globalOpts.dispSource
        %'GM'- Uses the grid method processing code developed by Grediac to
            %calculate displacements (the standard method)
        %'MatchID'- Optains displacements from match ID either grid method
            %or dic
quest='Choose source of displacement data';
btn1='Grid Images';
btn2='MatchID';
defbtn='Grid Images';
tempSource=questdlg(quest,'Displacement Source',...
    btn1,btn2,defbtn);

clear quest btn1 btn2 defbtn

switch tempSource
    case 'Grid Images'
        globalOpts.dispSource='GM';
        globalOpts.ExtractDICStrain=false;

    case 'MatchID'
        globalOpts.dispSource='MatchID';
        dispExt=questdlg('Use Match ID strains?');
        switch dispExt
            case 'Yes'
                globalOpts.ExtractDICStrain=true;
            case 'No'
                globalOpts.ExtractDICStrain=false;
            case 'Cancel'
                globalOpts.ExtractDICStrain=false;
        end

        clear dispExt
end

clear tempSource

%% Generate material properties structure
%Structure outputs
    %material-  with fields
        %name- Name of the material 
        %model- material model with options
            %Maxwell- Generalized maxwell modulus
        %rho- material density
        %Kinf- long term bulk modulus
        %Ginf- long term bulk modulus
        %nuinf- long term poisson's ratio
        %numEl- number of maxwell elements
        %rotAngle- (I'm not 100% sure what this field is)

prompt={'Material Name'
'Material Model'
'density (kg/m^3)'
'Kinf (Pa)'
    'Ginf (Pa)'
    'nuinf'
    'Number of Maxwell Elements'
    'rotation angle'};
%default inputs
defIn={'2El FE GM';'Maxwell';'1185';'2.0694e9';'1.1825e9';'0.26';'1';'0'};
dims=1;
tempProps=inputdlg(prompt,'Known Material Properties',dims,defIn);

material.name=char(cell2mat(tempProps(1)));
material.model=char(cell2mat(tempProps(2)));
material.rho=str2double(cell2mat(tempProps(3)));
material.Kinf=str2double(cell2mat(tempProps(4)));
material.Ginf=str2double(cell2mat(tempProps(5)));
material.nuinf=str2double(cell2mat(tempProps(6)));
material.numEl=str2double(cell2mat(tempProps(7)));
material.rotAngle=str2double(cell2mat(tempProps(8)));

clear prompt defIn dims tempProps

%% Generate SHear strain Smoothing opts
    %Structure output
        %smoothOpts.strain- with fields
            %spatialSmooth- True if spatially smoothing false if no
                %smoothing
            %spatialType- (Don't know exactly what this input is)
            %spatialAlgorithm- algorithm used to smooth for strain with
                %options
                    %'gaussian'- uses a guassian filter for the smoothing
            %spatialKernelSize- size of the smoothing window
            %spatialKernelStd- (I don't know what this is yet)
            %spatialEdgeMode-  (I don't know exactly what this is)
            %customEdgeMode-   (I dont know what this is)

Title='Strain Shear Smoothing Options';
prompt={'Shear Spatial Smoothing?';...
    'Type';... As of 2023/02/28 I need to know what this option means
     'Algorithm?';...
    'Kernel Size (px,px)';...
    'Kernel std';...
    'Edge Mode';...
    'Custome Edge Mode'};

defIn={'1';...
    'MATLAB';...
    'gaussian';...
    '11,11';...
    '10,10';...
    'replicate';...
    'extrapolate'};

dims=1;

tempParam=inputdlg(prompt,Title,dims,defIn);

Shear.smoothOpts.strain.spatialSmooth=str2double(cell2mat(tempParam(1)));
Shear.smoothOpts.strain.spatialSmooth=logical(...
    Shear.smoothOpts.strain.spatialSmooth);
Shear.smoothOpts.strain.spatialType=char(cell2mat(tempParam(2)));
Shear.smoothOpts.strain.spatialAlgorithm=char(cell2mat(tempParam(3)));
Shear.smoothOpts.strain.spatialKernelSize=str2num(cell2mat(tempParam(4)));
Shear.smoothOpts.strain.spatialKernelStd=str2num(cell2mat(tempParam(5)));
Shear.smoothOpts.strain.spatialEdgeMode=char(cell2mat(tempParam(6)));
Shear.smoothOpts.strain.customEdgeMode=char(cell2mat(tempParam(7)));

clear tempParam Title Prompt defIn dims

%% Bulk Smoothing Options 
Title='Strain Shear Smoothing Options';
prompt={'Bulk Spatial Smoothing?';...
    'Type';... As of 2023/02/28 I need to know what this option means
     'Algorithm?';...
    'Kernel Size (px,px)';...
    'Kernel std';...
    'Edge Mode';...
    'Custome Edge Mode'};

defIn={'1';...
    'MATLAB';...
    'gaussian';...
    '11,11';...
    '10,10';...
    'replicate';...
    'extrapolate'};

dims=1;

tempParam=inputdlg(prompt,Title,dims,defIn);

Bulk.smoothOpts.strain.spatialSmooth=str2double(cell2mat(tempParam(1)));
Bulk.smoothOpts.strain.spatialSmooth=logical(...
    Bulk.smoothOpts.strain.spatialSmooth);
Bulk.smoothOpts.strain.spatialType=char(cell2mat(tempParam(2)));
Bulk.smoothOpts.strain.spatialAlgorithm=char(cell2mat(tempParam(3)));
Bulk.smoothOpts.strain.spatialKernelSize=str2num(cell2mat(tempParam(4)));
Bulk.smoothOpts.strain.spatialKernelStd=str2num(cell2mat(tempParam(5)));
Bulk.smoothOpts.strain.spatialEdgeMode=char(cell2mat(tempParam(6)));
Bulk.smoothOpts.strain.customEdgeMode=char(cell2mat(tempParam(7)));

clear tempParam Title Prompt defIn dims
%% Generate Acceleratation smoothing options for Shear Identification
%Structure output
    %smoothOpts.acccel- with fields
        %temporalSmooth- true if smoothing false if not
        %temporalType- (I do not know what this means)
        %temporalAlgorithm- algorithm used to temporally smooth data
        %temporalKernelSize- size of the temporal smoothing kernal in
            %pixels
        %temporalKernelOrder- order of the svitsky golay filter

Title='Acceleration Smoothng Options for Shear';
prompt={'Temporal smoothing (1 or 0)';...
    'Type';...
    'algorithm';...
    'kernal size (frames,frames)';...
    'Filter order'};
defIn={'1';...
    'MATLAB';...
    'sgolay';...
    '11,11';...
    '3,3'};
dims=1;

tempParam=inputdlg(prompt,Title,dims,defIn);

Shear.smoothOpts.accel.temporalSmooth=str2double(cell2mat(tempParam(1)));
Shear.smoothOpts.accel.temporalSmooth=logical(...
    Shear.smoothOpts.accel.temporalSmooth);
Shear.smoothOpts.accel.temporalType=char(cell2mat(tempParam(2)));
Shear.smoothOpts.accel.temporalAlgorithm=char(cell2mat(tempParam(3)));
Shear.smoothOpts.accel.temporalKernelSize=str2num(cell2mat(tempParam(4)));
Shear.smoothOpts.accel.temporalKernelOrder=str2num(cell2mat(tempParam(5)));

clear Title prompt defIn dims tempParam

%% Generate Acceleratation smoothing options for Bulk Identification
%Structure output
    %smoothOpts.acccel- with fields
        %temporalSmooth- true if smoothing false if not
        %temporalType- (I do not know what this means)
        %temporalAlgorithm- algorithm used to temporally smooth data
        %temporalKernelSize- size of the temporal smoothing kernal in
            %pixels
        %temporalKernelOrder- order of the svitsky golay filter

Title='Acceleration Smoothng Options for Bulk';
prompt={'Temporal smoothing (1 or 0)';...
    'Type';...
    'algorithm';...
    'kernal size (frames,frames)';...
    'Filter order'};
defIn={'1';...
    'MATLAB';...
    'sgolay';...
    '11,11';...
    '3,3'};
dims=1;

tempParam=inputdlg(prompt,Title,dims,defIn);

Bulk.smoothOpts.accel.temporalSmooth=str2double(cell2mat(tempParam(1)));
Bulk.smoothOpts.accel.temporalSmooth=logical(...
    Bulk.smoothOpts.accel.temporalSmooth);
Bulk.smoothOpts.accel.temporalType=char(cell2mat(tempParam(2)));
Bulk.smoothOpts.accel.temporalAlgorithm=char(cell2mat(tempParam(3)));
Bulk.smoothOpts.accel.temporalKernelSize=str2num(cell2mat(tempParam(4)));
Bulk.smoothOpts.accel.temporalKernelOrder=str2num(cell2mat(tempParam(5)));

clear Title prompt defIn dims tempParam

%% Generate Shear Strain Rate Smoothing Options Structure
%structure Ouptut
    %smoothOpts.strainRate

Title='Shear Strain Rate Smoothing options';
prompt={'Temporal smoothing (1 or 0)';...
    'Type';...
    'algorithm';...
    'kernal size (px,px,px)';...
    'Filter order'};
defIn={'0';...
    'MATLAB';...
    'sgolay';...
    '11,11,11';...
    '3,3,3'};
dims=1;

tempParam=inputdlg(prompt,Title,dims,defIn);

Shear.smoothOpts.strainRate.temporalSmooth=str2double(cell2mat(tempParam(1)));
Shear.smoothOpts.strainRate.temporalSmooth=logical(...
    Shear.smoothOpts.strainRate.temporalSmooth);
Shear.smoothOpts.strainRate.temporalType=char(cell2mat(tempParam(2)));
Shear.smoothOpts.strainRate.temporalAlgorithm=char(cell2mat(tempParam(3)));
Shear.smoothOpts.strainRate.temporalKernelSize=str2num(cell2mat(tempParam(4)));
Shear.smoothOpts.strainRate.temporalKernelOrder=str2num(cell2mat(tempParam(5)));

clear Title prompt defIn dims tempParam

%% Generate Bulk Strain Rate Smoothing Options Structure
%structure Ouptut
    %smoothOpts.strainRate

Title='Bulk Strain Rate Smoothing options';
prompt={'Temporal smoothing (1 or 0)';...
    'Type';...
    'algorithm';...
    'kernal size (px,px,px)';...
    'Filter order'};
defIn={'0';...
    'MATLAB';...
    'sgolay';...
    '11,11,11';...
    '3,3,3'};
dims=1;

tempParam=inputdlg(prompt,Title,dims,defIn);

Bulk.smoothOpts.strainRate.temporalSmooth=str2double(cell2mat(tempParam(1)));
Bulk.smoothOpts.strainRate.temporalSmooth=logical(...
    Bulk.smoothOpts.strainRate.temporalSmooth);
Bulk.smoothOpts.strainRate.temporalType=char(cell2mat(tempParam(2)));
Bulk.smoothOpts.strainRate.temporalAlgorithm=char(cell2mat(tempParam(3)));
Bulk.smoothOpts.strainRate.temporalKernelSize=str2num(cell2mat(tempParam(4)));
Bulk.smoothOpts.strainRate.temporalKernelOrder=str2num(cell2mat(tempParam(5)));

clear Title prompt defIn dims tempParam


%% Set differentiation options
 %output structure
    %diffOpts.method- chooses method of temporal differentiation with
        %options
            %gradient- uses inbuilt matlab function gradient
            %cDiff- uses a custom central difference method

 quest='temporal differentiation method';
 Title='Temporal Differentiation Method';
 btn1='gradient';
 btn2='Central Difference';
 defbtn='gradient';
 
 tempOpts=questdlg(quest,Title,btn1,btn2,defbtn);

 switch tempOpts
     case 'gradient'
         diffOpts.method='gradient';
     case 'Central Difference'
         diffOpts.method='cDiff';
 end

 clear quest Title btn1 btn2 debtn tempOpts defbtn

 %% set displacement extrapolation options for Shear
 %output structure
    %extrapOpts.disp- Parameters for displacement corrections at the grid
        %edges using func_cropAndExtrapDispFields_v4 with fields:
            %fieldOpt- determines how to crop and extrapolate the
                %displacement fields
                    %XY- only crops the pixels that are directly lost (x on
                        %left and right edges, y on top and bottom)
                    %both- crops and extrapolates x and y displacements on
                        %all edges
            %cropPx1st- number of pixels to crop on the edges where the
                %grid ends
            %cropPx2nd- number of pixels to crop on  the edges where the
                %grid continues
            %extrapPx1st- number of pixels to correct on the edges where
                %the grid ends
            %extrapPx2nd- number of pixels to correct on the edges where
                %the gridcontinues
            %extrapMethod- order of the edge correction function with
                %options
                    %linear
                    %quadratic
            %extrapFitWinPx- number of pixels to use to obtain the best fit
                %for extrapolation
            %fixNaNKernel- number of pixels to fix NAN results in

Title='Displacement Extrapolation Options for Shear';
prompt={'Field option';...
    'Pixels to Crop on lost edges';...
    'Pixels to Crop on continuing edges';...
    'Pixels to Correct on lost edges';...
    'Pixels to Correct on continuing edges';...
    'Correction order';...
    'number of pixels to fit';...
    'Number of pixels to average to fix NaNs'};

defIn={'both';...
    '7';...
    '7';...
    '';...
    '';...
    'quadratic';...
    '12';...
    '18'};
dims=1;

tempOpts=inputdlg(prompt,Title,dims,defIn);

Shear.extrapOpts.disp.fieldOpt=char(cell2mat(tempOpts(1)));
Shear.extrapOpts.disp.cropPx1st=str2double(cell2mat(tempOpts(2)));
Shear.extrapOpts.disp.cropPx2nd=str2double(cell2mat(tempOpts(3)));
Shear.extrapOpts.disp.extrapPx1st=str2double(cell2mat(tempOpts(4)));
Shear.extrapOpts.disp.extrapPx2nd=str2double(cell2mat(tempOpts(5)));
Shear.extrapOpts.disp.extrapMethod=char(cell2mat(tempOpts(6)));
Shear.extrapOpts.disp.extrapFitWinPx=str2double(cell2mat(tempOpts(7)));
Shear.extrapOpts.disp.fixNaNKernel=str2double(cell2mat(tempOpts(8)));

clear Title prompt defIn dims tempOpts

 %% set displacement extrapolation options for Bulk
 %output structure
    %extrapOpts.disp- Parameters for displacement corrections at the grid
        %edges using func_cropAndExtrapDispFields_v4 with fields:
            %fieldOpt- determines how to crop and extrapolate the
                %displacement fields
                    %XY- only crops the pixels that are directly lost (x on
                        %left and right edges, y on top and bottom)
                    %both- crops and extrapolates x and y displacements on
                        %all edges
            %cropPx1st- number of pixels to crop on the edges where the
                %grid ends
            %cropPx2nd- number of pixels to crop on  the edges where the
                %grid continues
            %extrapPx1st- number of pixels to correct on the edges where
                %the grid ends
            %extrapPx2nd- number of pixels to correct on the edges where
                %the gridcontinues
            %extrapMethod- order of the edge correction function with
                %options
                    %linear
                    %quadratic
            %extrapFitWinPx- number of pixels to use to obtain the best fit
                %for extrapolation
            %fixNaNKernel- number of pixels to fix NAN results in

Title='Displacement Extrapolation Options for Bulk';
prompt={'Field option';...
    'Pixels to Crop on lost edges';...
    'Pixels to Crop on continuing edges';...
    'Pixels to Correct on lost edges';...
    'Pixels to Correct on continuing edges';...
    'Correction order';...
    'number of pixels to fit';...
    'Number of pixels to average to fix NaNs'};

defIn={'both';...
    '7';...
    '7';...
    '';...
    '';...
    'quadratic';...
    '12';...
    '18'};
dims=1;

tempOpts=inputdlg(prompt,Title,dims,defIn);

Bulk.extrapOpts.disp.fieldOpt=char(cell2mat(tempOpts(1)));
Bulk.extrapOpts.disp.cropPx1st=str2double(cell2mat(tempOpts(2)));
Bulk.extrapOpts.disp.cropPx2nd=str2double(cell2mat(tempOpts(3)));
Bulk.extrapOpts.disp.extrapPx1st=str2double(cell2mat(tempOpts(4)));
Bulk.extrapOpts.disp.extrapPx2nd=str2double(cell2mat(tempOpts(5)));
Bulk.extrapOpts.disp.extrapMethod=char(cell2mat(tempOpts(6)));
Bulk.extrapOpts.disp.extrapFitWinPx=str2double(cell2mat(tempOpts(7)));
Bulk.extrapOpts.disp.fixNaNKernel=str2double(cell2mat(tempOpts(8)));

clear Title prompt defIn dims tempOpts

%% Set Shear strain extrapolation options 
 %output structure
    %extrapOpts.strain- Parameters for strain cropping and corrections at 
        %the grid edges using func_cropStrainFields and 
        %func_extrapolateStrainFields with fields:
            %postExtrapOn- logical value to determine whether strains are
                %extrapolated after displacement corrections
                    %true- strain is cropped and extrapolated
                    %false- strain only uses corrected displacement values
            %fieldOpt- determines how to crop and extrapolate the
                %displacement fields
                    %XY- only crops the pixels that are directly lost (x on
                        %left and right edges, y on top and bottom)
                    %both- crops and extrapolates x and y displacements on
                        %all edges
            %extrapPx1st- number of pixels to crop and extrapolate on the 
                %edges where the grid ends
            %extrapPx2nd- number of pixels to crop and extrapolate on the
                %edges where the grid continues
            %extrapMethod- order of the edge correction function with
                %options
                    %linear
                    %quadratic
            %extrapFitWinPx- number of pixels to use to obtain the best fit
                %for extrapolation
            %enforceGlobBCs- determines which strain boundary conditions to
                %enforce
                    %'shear'- extrapolates shear strain to 0 on the edges
                    %no other options currently
            %enforecGlobBCsPx- number of pixels to enforce the shear
                %boundary condition on

Title='Strain Extrapolation Options for Shear';
prompt={'Correct strain?'
    'Field option';...
    'Pixels to Correct on lost edges';...
    'Pixels to Correct on continuing edges';...
    'Correction order';...
    'number of pixels to fit';...
    'Which global boundary conditions to fix';...
    'Number of pixels enforce shear BC'};

defIn={'1'
    'both';...
    '7';...
    '7';...
    'quadratic';...
    '12';...
    'shear'
    '7'};
dims=1;

tempOpts=inputdlg(prompt,Title,dims,defIn);


Shear.extrapOpts.strain.postExtrapOn=str2double(cell2mat(tempOpts(1)));
Shear.extrapOpts.strain.postExtrapOn=logical(...
    Shear.extrapOpts.strain.postExtrapOn);
Shear.extrapOpts.strain.fieldOpt=char(cell2mat(tempOpts(2)));
Shear.extrapOpts.strain.extrapPx1st=str2double(cell2mat(tempOpts(3)));
Shear.extrapOpts.strain.extrapPx2nd=str2double(cell2mat(tempOpts(4)));
Shear.extrapOpts.strain.extrapMethod=char(cell2mat(tempOpts(5)));
Shear.extrapOpts.strain.extrapFitWinPx=str2double(cell2mat(tempOpts(6)));
Shear.extrapOpts.strain.enforceGlobBCs=char(cell2mat(tempOpts(7)));
Shear.extrapOpts.strain.enforceGlobBCsPx=str2double(cell2mat(tempOpts(8)));

clear Title prompt defIn dims tempOpts

%% Set Bulk strain extrapolation options 
 %output structure
    %extrapOpts.strain- Parameters for strain cropping and corrections at 
        %the grid edges using func_cropStrainFields and 
        %func_extrapolateStrainFields with fields:
            %postExtrapOn- logical value to determine whether strains are
                %extrapolated after displacement corrections
                    %true- strain is cropped and extrapolated
                    %false- strain only uses corrected displacement values
            %fieldOpt- determines how to crop and extrapolate the
                %displacement fields
                    %XY- only crops the pixels that are directly lost (x on
                        %left and right edges, y on top and bottom)
                    %both- crops and extrapolates x and y displacements on
                        %all edges
            %extrapPx1st- number of pixels to crop and extrapolate on the 
                %edges where the grid ends
            %extrapPx2nd- number of pixels to crop and extrapolate on the
                %edges where the grid continues
            %extrapMethod- order of the edge correction function with
                %options
                    %linear
                    %quadratic
            %extrapFitWinPx- number of pixels to use to obtain the best fit
                %for extrapolation
            %enforceGlobBCs- determines which strain boundary conditions to
                %enforce
                    %'shear'- extrapolates shear strain to 0 on the edges
                    %no other options currently
            %enforecGlobBCsPx- number of pixels to enforce the shear
                %boundary condition on

Title='Strain Extrapolation Options for Bulk';
prompt={'Correct strain?'
    'Field option';...
    'Pixels to Correct on lost edges';...
    'Pixels to Correct on continuing edges';...
    'Correction order';...
    'number of pixels to fit';...
    'Which global boundary conditions to fix';...
    'Number of pixels enforce shear BC'};

defIn={'1'
    'both';...
    '7';...
    '7';...
    'quadratic';...
    '12';...
    'shear'
    '7'};
dims=1;

tempOpts=inputdlg(prompt,Title,dims,defIn);


Bulk.extrapOpts.strain.postExtrapOn=str2double(cell2mat(tempOpts(1)));
Bulk.extrapOpts.strain.postExtrapOn=logical(Bulk.extrapOpts.strain.postExtrapOn);
Bulk.extrapOpts.strain.fieldOpt=char(cell2mat(tempOpts(2)));
Bulk.extrapOpts.strain.extrapPx1st=str2double(cell2mat(tempOpts(3)));
Bulk.extrapOpts.strain.extrapPx2nd=str2double(cell2mat(tempOpts(4)));
Bulk.extrapOpts.strain.extrapMethod=char(cell2mat(tempOpts(5)));
Bulk.extrapOpts.strain.extrapFitWinPx=str2double(cell2mat(tempOpts(6)));
Bulk.extrapOpts.strain.enforceGlobBCs=char(cell2mat(tempOpts(7)));
Bulk.extrapOpts.strain.enforceGlobBCsPx=str2double(cell2mat(tempOpts(8)));

clear Title prompt defIn dims tempOpts

%% Set acceleration extrapoloation options for Shear
%output structure
    %extrapOpts.accel Parameters for temporal padding of displacements 
    %using func_temporalPadDisp along with acceleration cropping and 
    %correctionsa t the grid edges using func_cropAccelFields and 
        %func_extrapolateAccelFields with fields:
            %tempPadOn- logical value to determine if temporal padding is
                %performed. Moved from smoothingOpts from before.
                    %true- temporal padding is performed
                    %false- no temporal padding
            %tempPadFrames- vector specifying number of frames to pad disp
                %[# of frames to pad at beginning, # of frames at end]
            %tempPadMethod- Order of temporal padding
                %replicate- default option. Just coppies first and end
                    %frames
                %quadratic- quadratically fits disp fields temporally
                %linear- linearly fits and extrapolates disp temporally
            %tempPadFitWin- vector of fitting window for temporal padding
                %[# frames to fit at start, # frames to fit at end]
            %postSpatExtrapOn- logical value to determine whether
                %accelerations are extrapolated on the edges after disp 
                    %true- accel is cropped and extrapolated
                    %false- accel only uses corrected and smoothed disp
            %fieldOpt- determines how to crop and extrapolate the
                %displacement fields
                    %XY- only crops the pixels that are directly lost (x on
                        %left and right edges, y on top and bottom)
                    %both- crops and extrapolates x and y displacements on
                        %all edges
            %extrapPx1st- number of pixels to crop and extrapolate on the 
                %edges where the grid ends
            %extrapPx2nd- number of pixels to crop and extrapolate on the
                %edges where the grid continues
            %extrapMethod- order of the edge correction function with
                %options
                    %linear
                    %quadratic
            %extrapFitWinPx- number of pixels to use to obtain the best fit
                %for extrapolation

Title='Acceleration Extrapolation and padding Options for Shear';
prompt={'Pad Temporally?';...
    '# frames to pad (beginning,end)'
    'Padding Method'
    'Temporal fit window (beginning,end)'
    'Correct accelerations? 1 or 0';...
    'Field option';...
    'Pixels to Correct on lost edges';...
    'Pixels to Correct on continuing edges';...
    'Correction order';...
    'number of pixels to fit'};

defIn={'1'
    '5,0';...
    'replicate'
    '11,11'
    '1'
    'both'
    '7';...
    '7';...
    'quadratic';...
    '12'};
dims=1;

tempOpts=inputdlg(prompt,Title,dims,defIn);

Shear.extrapOpts.accel.tempPadOn=str2double(cell2mat(tempOpts(1)));
Shear.extrapOpts.accel.tempPadOn=logical(Shear.extrapOpts.accel.tempPadOn);
Shear.extrapOpts.accel.tempPadFrames=str2num(cell2mat(tempOpts(2)));
Shear.extrapOpts.accel.tempPadMethod=char(cell2mat(tempOpts(3)));
Shear.extrapOpts.accel.tempPadFitWin=str2num(cell2mat(tempOpts(4)));
Shear.extrapOpts.accel.postSpatExtrapOn=str2double(cell2mat(tempOpts(5)));
Shear.extrapOpts.accel.postSpatExtrapOn=logical(...
    Shear.extrapOpts.accel.postSpatExtrapOn);
Shear.extrapOpts.accel.fieldOpt=char(cell2mat(tempOpts(6)));
Shear.extrapOpts.accel.extrapPx1st=str2double(cell2mat(tempOpts(7)));
Shear.extrapOpts.accel.extrapPx2nd=str2double(cell2mat(tempOpts(8)));
Shear.extrapOpts.accel.extrapMethod=char(cell2mat(tempOpts(9)));
Shear.extrapOpts.accel.extrapFitWinPx=str2double(cell2mat(tempOpts(10)));

clear Title prompt defIn dims tempOpts

%% Set acceleration extrapoloation options for Bulk Identification
%output structure
    %extrapOpts.accel Parameters for temporal padding of displacements 
    %using func_temporalPadDisp along with acceleration cropping and 
    %correctionsa t the grid edges using func_cropAccelFields and 
        %func_extrapolateAccelFields with fields:
            %tempPadOn- logical value to determine if temporal padding is
                %performed. Moved from smoothingOpts from before.
                    %true- temporal padding is performed
                    %false- no temporal padding
            %tempPadFrames- vector specifying number of frames to pad disp
                %[# of frames to pad at beginning, # of frames at end]
            %tempPadMethod- Order of temporal padding
                %replicate- default option. Just coppies first and end
                    %frames
                %quadratic- quadratically fits disp fields temporally
                %linear- linearly fits and extrapolates disp temporally
            %tempPadFitWin- vector of fitting window for temporal padding
                %[# frames to fit at start, # frames to fit at end]
            %postSpatExtrapOn- logical value to determine whether
                %accelerations are extrapolated on the edges after disp 
                    %true- accel is cropped and extrapolated
                    %false- accel only uses corrected and smoothed disp
            %fieldOpt- determines how to crop and extrapolate the
                %displacement fields
                    %XY- only crops the pixels that are directly lost (x on
                        %left and right edges, y on top and bottom)
                    %both- crops and extrapolates x and y displacements on
                        %all edges
            %extrapPx1st- number of pixels to crop and extrapolate on the 
                %edges where the grid ends
            %extrapPx2nd- number of pixels to crop and extrapolate on the
                %edges where the grid continues
            %extrapMethod- order of the edge correction function with
                %options
                    %linear
                    %quadratic
            %extrapFitWinPx- number of pixels to use to obtain the best fit
                %for extrapolation

Title='Acceleration Extrapolation and padding Options for Bulk';
prompt={'Pad Temporally?';...
    '# frames to pad (beginning,end)'
    'Padding Method'
    'Temporal fit window (beginning,end)'
    'Correct accelerations? 1 or 0';...
    'Field option';...
    'Pixels to Correct on lost edges';...
    'Pixels to Correct on continuing edges';...
    'Correction order';...
    'number of pixels to fit'};

defIn={'1'
    '5,0';...
    'replicate'
    '11,11'
    '1'
    'both'
    '7';...
    '7';...
    'quadratic';...
    '12'};
dims=1;

tempOpts=inputdlg(prompt,Title,dims,defIn);

Bulk.extrapOpts.accel.tempPadOn=str2double(cell2mat(tempOpts(1)));
Bulk.extrapOpts.accel.tempPadOn=logical(Bulk.extrapOpts.accel.tempPadOn);
Bulk.extrapOpts.accel.tempPadFrames=str2num(cell2mat(tempOpts(2)));
Bulk.extrapOpts.accel.tempPadMethod=char(cell2mat(tempOpts(3)));
Bulk.extrapOpts.accel.tempPadFitWin=str2num(cell2mat(tempOpts(4)));
Bulk.extrapOpts.accel.postSpatExtrapOn=str2double(cell2mat(tempOpts(5)));
Bulk.extrapOpts.accel.postSpatExtrapOn=logical(...
    Bulk.extrapOpts.accel.postSpatExtrapOn);
Bulk.extrapOpts.accel.fieldOpt=char(cell2mat(tempOpts(6)));
Bulk.extrapOpts.accel.extrapPx1st=str2double(cell2mat(tempOpts(7)));
Bulk.extrapOpts.accel.extrapPx2nd=str2double(cell2mat(tempOpts(8)));
Bulk.extrapOpts.accel.extrapMethod=char(cell2mat(tempOpts(9)));
Bulk.extrapOpts.accel.extrapFitWinPx=str2double(cell2mat(tempOpts(10)));

clear Title prompt defIn dims tempOpts

%% Generate genSSCurveOpts
%Not needed for current coding tasks. To be implemented later so that I can
    %integrate my work with existing IBIII codes

%% Generate global options
%Note: Many of these parameters aren't used in the current iteration of the
    %Viscoelastic IBII codes. This section will be completed once
    %viscoelastic IBII is implemented with the existing IBII codes

%output structure
    %global opts- options governing the global parameter identification
        %options with fields
            %imageDefAddNoise- Add noise to deformed images, logical
            %smoothingOn- logical deterines if smoothing
            %matModel- material model
            %fieldComponents-
            %hardCodeFreeEdge- how free edge is determined logical
            %freeEdge- determines where the free edge is in the images
                %left- specimen is impacted on the right
                %right- specimen is impacted on the left
            %plotImages- determines whether to evaluation images
            %plotImageSeqFlags-
            %plotImageSeqMatFlags-
            %processOptVF- Run IBII with optimized virtual fields, logical
            %processShearSG- Whether or not to run the shear stress gage,
                %logical
            %processGenSSCurves- Use generalized stress--strain curves,
                %logical
            %identifyStrength- Whether or not to identify strength (if
                %specimen fails), logical
            %calcEBal- logical 
            %plotPresentationImages- whether to plot images for
                %presentations, logical
            %reduceMemory- whether to keep or delete intermediate data to
                %save on RAM usage
            %saveProcData- whether or not to save processed data, string
                %'yes' saves the data

           
Title='Global Processing Options';
prompt={'Add Noise to Image Deformation';...
    'Smooth?'
    'Material Model'
    'field components'
    'hardcode free edge';...
    'Free edge';...
    'Plot images';...
    'plot image flags'
    'plot image Mat flags'
    'process Optimized virtual fields';...
    'Process Shear Stress Gauge';...
    'Generalized Stress--strain'
    'identify strength'
    'calc E Balance'
    'plot presentation imnages'
    'Reduce Memory'
    'save Processed data'};

defIn={'0'
    '1'
    'Maxwell'
    'all'
    '1'
    'Left'
    'prompt'
    '1,1,0,1'
    '0,0,0'
    '0'
    '1'
    '0'
    '0'
    '0'
    '0'
    '1'
    'yes'};

dims=1;

tempOpts=inputdlg(prompt,Title,dims,defIn);

globalOpts.imageDefAddNoise=str2double(cell2mat(tempOpts(1)));
globalOpts.imageDefAddNoise=logical(globalOpts.imageDefAddNoise);
globalOpts.smoothingOn=str2double(cell2mat(tempOpts(2)));
globalOpts.smoothingOn=logical(globalOpts.smoothingOn);
globalOpts.matModel=char(cell2mat(tempOpts(3)));
globalOpts.fieldComponents=char(cell2mat(tempOpts(4)));
globalOpts.hardCodeFreeEdge=str2double(cell2mat(tempOpts(5)));
globalOpts.hardCodeFreeEdge=logical(globalOpts.hardCodeFreeEdge);
globalOpts.freeEdge=char(cell2mat(tempOpts(6)));
globalOpts.plotImages=char(cell2mat(tempOpts(7)));
globalOpts.plotImageSeqFlags=str2num(cell2mat(tempOpts(8)));
globalOpts.plotImageSeqMatFlags=str2num(cell2mat(tempOpts(9)));
globalOpts.processOptVF=str2double(cell2mat(tempOpts(10)));
globalOpts.processOptVF=logical(globalOpts.processOptVF);
globalOpts.processShearSG=str2double(cell2mat(tempOpts(11)));
globalOpts.processShearSG=logical(globalOpts.processShearSG);
globalOpts.processGenSSCurves=str2double(cell2mat(tempOpts(12)));
globalOpts.processGenSSCurves=logical(globalOpts.processGenSSCurves);
globalOpts.identifyStrength=str2double(cell2mat(tempOpts(13)));
globalOpts.identifyStrength=logical(globalOpts.identifyStrength);
globalOpts.calcEBal=str2double(cell2mat(tempOpts(14)));
globalOpts.calcEBal=logical(globalOpts.calcEBal);
globalOpts.plotPresentationImages=str2double(cell2mat(tempOpts(15)));
globalOpts.plotPresentationImages=...
    logical(globalOpts.plotPresentationImages);
globalOpts.reduceMemory=str2double(cell2mat(tempOpts(16)));
globalOpts.reduceMemory=logical(globalOpts.reduceMemory);
globalOpts.saveProcData=char(cell2mat(tempOpts(17)));

clear Title prompt defIn dims tempOpts

%% Input Grid Parameters
switch globalOpts.dispSource
    case 'GM'
        %output structure
        %grid- structure containing grid parameters with fields
        %name-string containing the name of the particular grid
        %pitch- grid pitch (m)
        %pxPerPeriod- number of pixels per grid period
        %rotAngle-
        %length- length of the specimen in ROI m
        %height- height of the specimen ROI in m
        %mPerPx- length of a pixel in m
        %numXPeriods- number of periods along the length of the specimen
        %asymmPitch- is the gid asymmetrical

        Title='Grid Parameters';
        prompt={'name'
            'pitch'
            '# pixels sampling a period'
            'rotation angle'
            'Specimen Length (IN FOV)'
            'specimen height (IN FOV)'
            'asymmetric Pitch?'};
        defIn={'09mm_5pxPP_printed_RapidPress'
            '0.9e-3'
            '5'
            '0'
            '70e-3'
            '44e-3'
            '0'};
        dims=1;
        tempParam=inputdlg(prompt,Title,dims,defIn);

        grid.name=char(cell2mat(tempParam(1)));
        grid.pitch=str2double(cell2mat(tempParam(2)));
        grid.pxPerPeriod=str2double(cell2mat(tempParam(3)));
        grid.rotAngle=str2double(cell2mat(tempParam(4)));
        grid.length=str2double(cell2mat(tempParam(5)));
        grid.height=str2double(cell2mat(tempParam(6)));
        grid.asymmPitch=str2double(cell2mat(tempParam(7)));

        %Calculate remaining grid parameters
        grid.mPerPx=grid.pitch/grid.pxPerPeriod;
        grid.numXPeriods=grid.length/grid.pitch;

        
    case 'MatchID'
        prompt={'SS size (px)';'ST (px)';'DIC displacement units'};
        Title='DIC Parameters';
        defIn={'25','3','1e-3'};
        dims=1;
        opts.Resize='on';
        opts.WindowStyle='normal';
       tempParam=inputdlg(prompt,Title,dims,defIn,opts);
       DIC.ST=str2double(tempParam(2));
       DIC.SS=str2double(tempParam(1));
       DIC.ConvUnit=str2double(tempParam(3));
end

clear Title prompt defIn dims tempParam

%%  set grid method processing options
%Note I have never chaged these so they are hardcoded

gridMethodOpts.windowFlag=1;
gridMethodOpts.windowWidth=1;
gridMethodOpts.dispCalcMethod=2;
gridMethodOpts.temporalUnwrap=1;
gridMethodOpts.debug=0;
gridMethodOpts.hardCodeRotAngle=1;
gridMethodOpts.autoLoadProccessedDataFile=0;

%% Image Noise 
%output structure
    %imageNoise- parameters for noise added to the images with fiedls
        %addNoise- logical true for noise false for no noise
        %pcNoise- level of noise to ad (% of dynamic range)
        %bits- number of bits encoding the image
        %convtoUInt16- 

quest='Add Noise to Images';
NoiseChoice=questdlg(quest);
switch NoiseChoice
    case 'No'
        imageNoise.addNoise=false;
    case 'Yes'
        imageNoise.addNoise=true;
        Title='Noise Parameters';
        prompt={'noise level (% of dynamic range)'
            'bits'
            'convToUInt16'};
        defIn={'0.4'
            '16'
            'true'};
        tempOpts=inputdlg(prompt,Title,1,defIn);
       
        imageNoise.pcNoise=str2double(cell2mat(tempOpts(1)));
        imageNoise.bits=str2double(cell2mat(tempOpts(2)));
        imageNoise.convToUInt16=str2num(cell2mat(tempOpts(3)));

        clear  prompt Title defIn tempOpts
end 
clear NoiseChoice quest

%% Set Specimen Properties
%output structure
    %specimen- structure containing specimen properties
        %length- length of specimen (total length) (m) 
        %height- height of specimen (total height) (m)
        %thickness- specimen thickness in m
        %volume- volume of the specimen 
        %mass- mass of the specimen
        %freeEdge- specimen free edge

%Make sure free edge matches global free edge
specimen.freeEdge=globalOpts.freeEdge;

% input other specimen properties
Title='Specimen Properties';
prompt={'Total specimen length (m)'
    'Total specimen height (m)'
    'measured specimen thickness (m)'
    'measured specimen mass (kg)'};
defIn={'70e-3'
    '44e-3'
    '5e-3'
    '0.018249'};
tempProps=inputdlg(prompt,Title,1,defIn);

specimen.length=str2double(cell2mat(tempProps(1)));
specimen.height=str2double(cell2mat(tempProps(2)));
specimen.thickness=str2double(cell2mat(tempProps(3)));
specimen.volume=specimen.length*specimen.height*specimen.thickness;
specimen.mass=str2double(cell2mat(tempProps(4)));

%Update density to use measured properties
%material.rho=specimen.mass/specimen.volume;


clear Title prompt defIn tempProps

%% create time data structure
%output structure
    %time- time properties of the experiment with fields
        %frameRate- camera frame rate in fps
        %step- time between frames
        %cutFrames- 
        %numFrames- total number of recorded frames
        %vector of time increments
 
Title='Time Properties';
prompt={'Frame Rate (fps)'
    'Record Length (frames)'
    'Cut Frames'};
defIn={'2e6'
    '128'
    '1'};
tempProps=inputdlg(prompt,Title,1,defIn);

time.frameRate=str2double(cell2mat(tempProps(1)));
time.step=1/time.frameRate;
time.numFrames=str2double(cell2mat(tempProps(2)));
time.cutFrames=str2double(cell2mat(tempProps(3)));
time.vec=0:time.step:(time.numFrames-1)*time.step;

clear Title prompt defIn tempProps


%% Create CondOpts data structure to condition data before pasing into the
    %cost function

    Title='Data conditioning options';
    prompt={'Downsample temporally (1 or 0)?'
        'Temporal downsampling factor (frames)'
        'Frames to remove from the end'
        'Frames to remove from the beginning'
        'Downsampling factor in the X'
        'Downsampling factor in the Y'
        'Pixels to remove on impact edge'
        'Pixels to remove on Free Edge'};
   defIn={'0'
       '1'
       '4'
       '0'
       '3'
       '1'
       '10'
       '50'};
  dlgOpts.Resize='on';
  dlgOpts.WindowStyle='normal';
 
  tempOpts=inputdlg(prompt,Title,1,defIn,dlgOpts);

  CondOpts.TempDS=logical(str2double(cell2mat(tempOpts(1))));
  CondOpts.Tds=str2double(cell2mat(tempOpts(2)));
  CondOpts.CutEndFrames=str2double(cell2mat(tempOpts(3)));
  CondOpts.CutStartFrames=str2double(cell2mat(tempOpts(4)));
  CondOpts.Xds=str2double(cell2mat(tempOpts(5)));
  CondOpts.Yds=str2double(cell2mat(tempOpts(6)));
  CondOpts.ImpCens=str2double(cell2mat(tempOpts(7)));
  CondOpts.FreeCens=str2double(cell2mat(tempOpts(8)));


  clear tempOpts Title prompt defIn



%% Generate SolveOpts structure

Title='Cost function solving options';
prompt={'Weight of the axial cost function (fraction of 1)'
    'Constant Poisson Ratio? (0 or 1)'
    'KGsame (0 or 1)' %this notation is confusing may remove in future iterations
    'identify long term moduli (1 or 0)'
    'which moduli to identify'
    'minimization function'};
defIn={'0.5'
    '0'
    '0'
    '0'
    'KG'
    'fimincon'};
  dlgOpts.Resize='on';
  dlgOpts.WindowStyle='normal';
 
  tempOpts=inputdlg(prompt,Title,1,defIn,diffOpts);

  SolveOpts.wK=str2double(cell2mat(tempOpts(1)));
  SolveOpts.wG=1-SolveOpts.wK;
  SolveOpts.constnu=str2double(cell2mat(tempOpts(2)));
  SolveOpts.KGsame=str2double(cell2mat(tempOpts(3)));
  SolveOpts.identEinf=str2double(cell2mat(tempOpts(4)));
  SolveOpts.identForm=char(cell2mat(tempOpts(5)));
  SolveOpts.minFunc=char(cell2mat(tempOpts(6)));


  clear prompt Title defIn tempOpts


%% minOpts data structure
Title='Solver options';
prompt={'Algorithm'
    'Compute in parallel? (0 or 1)'
    'step tolerance'};
defIn={'interior-point'
    '1'
    '1e-3'};
dlgOpts.Resize='on';
dlgOpts.WindowStyle='normal';

tempOpts=inputdlg(prompt,Title,1,defIn,dlgOpts);

minOpts.Algorithm=char(cell2mat(tempOpts(1)));
minOpts.UseParallel=str2double(cell2mat(tempOpts(2)));
minOpts.StepTolerance=str2double(cell2mat(tempOpts(3)));

clear prompt Title defIn tempOpts
%% Stuff to add
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Structures to implement
    %stressGaugeOpts
    %stressGaugeShearOpts
    %VFOpts
    %VFPlotOpts



%% Save Processing Parameters
saveName=strcat(SaveDir,'/',procDeg,'_ProcessParameters_V4.mat');
clear SaveDir
save(saveName)

