function [OutputStruct] = func_subGMEdgeFields(time,FEdisp,...
GMdisp,GMstrain,GMaccel,pos,grid,...
    saveDir,Desig,...
    subOpts,smoothingOpts,globalOpts,extrapOpts,diffOpts)

%This function is written to substitute interpolated finite element shear
    %on to two grid pitches from the top and bottom free surfaces of an
    %IBII specimen. Finite element displacements interpolated onto the grid
    %coordinate system are input along with the originally extracted grid
    %method data. The interpolated finite element data is then substituted
    %onto the 4 specimen edges according to the method specified in the
    %inputs, and the strain and acceleration fields are recalculated. This
    %script outputs a structure containing the original inputs along with
    %the substituted displacement and recalulated accelerations and
    %strains. The function also saves the data in a separate file to be
    %reffered to later.

%Author: Andrew Matejunas
%Date Created: 2022-09-26

%Change Log:
    %2022-09-26: Initial version created

%Function Inputs:
    %time- structure containing time information
    %FEdisp- Displacements obtained from pure finite elemenents
        %interpolated onto grid method coordinates
    %GMdisp- Displacements obtained from the grid method processing of
        %synthetic grid images with edge errors
    %GMstrain- orignially calculated strain with edge errors
    %GMaccel- originally calculated accelerations with edge errors
    %pos-  Structure Containing x and y coordinates from grid mefthod
        %data
    %subOpts- structure containing substitution options with fields
        %method- Method used to determine what pixels to substitute
            %'GridPeriod'- substitution is performed on a certain number of
                %grid pitches
            %'SpatialKernal- Substitution is performed over a
                %spatial smoothing kernal
        %PitchNum- if subOpts.method='GridPeriod' specifies number of grid
            %pitches over which substitution will be perfomed
   %smoothingOpts- standard structure containing instructions on how to
        %smooth the displacement fields spatially and temporally before
        %strain and acceleration calucalation
   %globalOpts- global options used in the strain calculations
   %extrapOps- extrapolation options used in the 
   %diffOpts- temporal differentiation options used in accccelration

        

%Function Outputs:
    %OutputStruc- structure containing output kinematic fields with fields
       %time- Records time structure
       %FEdisp- Records orignial FE displacement structure
       %GMdisp- Records initial grid method displacement with errors
       %GMstrain- records intial grid method strains with errors
       %GMaccel- records initial grid method accelerations with errors
       %disp- output displacement fields with edge inteterpolation
       %accel- outputs acceleration fields with edge interpolation
       %strain- outputs strain fields with edge interpolation
       %subOpts- Options used to perform substitution
       %smoothingOpts- record of smoothing parameters
       %globalOpts- global options used in the strain calculations
       %extrapOps- extrapolation options used in the 
       %diffOpts- temporal differentiation options used in accccelration
       %SaveFile- file containing the original and substituted fields 
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Substitute edge pixels
switch subOpts.method
    case 'GridPeriod'
        %define interval of pixels over the substitution will be performed
        SubIntY=subOpts.PitchNum*grid.pxPerPeriod;
        SubIntX=SubIntY;
   
   case 'SpatialKernal'
        if smoothingOpts.spatialKernal(1)>=grid.pxPerPeriod
           SubIntX=smoothingOpts.spatialKernal(1);
        else 
            SubIntX=grid.pxPerPeriod;
        end

        if smoothingOpts.spatialkernal(2)>=grid.pxPerPeriod
            SubIntY=smoothingOpts.spatialKernal;
        else
            SubIntY=grid.pxPerPeriod;
        end
end

%% Define indexes to substitute
 NumY=length(pos.y);
 NumX=length(pos.x);

 %top edge
 subInTop=(NumY-SubIntY):NumY;
        
 %bottom edge
 subInBot=1:SubIntY+1;

 %left edge
 subInLeft=1:SubIntX+1;

 %Right edge
 subInRight=(NumX-SubIntX):NumX;

%% Perform substitution
fprintf('Performing Edge substitutions \n')
%start with extracted grid method displacements
disp.x=GMdisp.x;
disp.y=GMdisp.y;

%Substitute top Edge
disp.x(subInTop,:,:)=FEdisp.x(subInTop,:,:);
disp.y(subInTop,:,:)=FEdisp.y(subInTop,:,:);
%bottom edge
disp.x(subInBot,:,:)=FEdisp.x(subInBot,:,:);
disp.y(subInBo,:,:)=FEdisp.y(subInBot,:,:);
%left edge
disp.x(:,subInLeft,:)=FEdisp.x(:,subInLeft,:);
disp.y(:,subInLeft,:)=FEdisp.y(:,subInLeft,:);
%right edge
disp.x(:,subInRight,:)=FEdisp.x(:,subInRight,:);
disp.y(:,subInRight,:)=FEdisp.y(:,subInRight,:);

%% Calculate Strains from displacements
fprintf('Calculating strains from substituted displacement fields \n')
[disp,strain,strainRate] = func_smoothCalcStrain(globalOpts,pos,time,...
            grid,disp,smoothingOpts,extrapOpts);
%% Calculate Accelerations from substituted displacements
fprintf('Calculating acceleration from siustituted displacement fields.\n')
 [disp,~,accel] = func_smoothCalcAccel(pos,time,grid,disp,smoothingOpts,...
            extrapOpts,diffOpts);

 %% Save Data
 fprintf('Saving substituted data \n')
 SaveFile=strcat(saveDir,'/',Desig,'_',subOpts.method, ...
     '_EdgeSubStitution.mat');
 save(SaveFile)
 OutputStruct.SaveFile=SaveFile;

%% ensure inputs are also output
OutputStruct.FEdisp=FEdisp;
OutputStruct.GMdisp=GMdisp;
OutputStruct.GMstrain=GMstrain;
OutputStruct.GMaccel=GMaccel;
OutputStruct.time=time;
OutputStruct.pos=pos;
OutputStruct.grid=grid;
OutputStruct.subOpts=subOpts;
OutputStruct.smoothingOpts=smoothingOpts;
OutputStruct.globalOpts=globalOpts;
OutputStruct.extrapOpts=extrapOpts;
OutputStruct.diffOpts=diffOpts;

%% Output the displacement, strain, strain rate, and acceleration fields
OutputStruct.disp=disp;
OutputStruct.strain=strain;
OutputStruct.strainRate=strainRate;
OutputStruct.accel=accel;

end