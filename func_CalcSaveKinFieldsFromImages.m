function [OutputStruct] =...
    func_CalcSaveKinFieldsFromImagesNoSmooth(MainSaveDir,...
    imageFile,ImagePath,ParentDesig,NoiseCall, ...
        grid,gridMethodOpts,imageNoise,specimen,DispCorr,globalOpts)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%% Set Smoothing options
%Hardcode No smoothing
smoothingOpts.spatialSmooth=0;
smoothingOpts.WATempSmooth=0;
smoothingOpts.FFTempSmooth=0;


%% Run grid method processing algorithm without noise
[grid,pos,disp] = func_gridMethodImageProcessing_AJM(imagePath,...
                imageFile,...
                grid,gridMethodOpts,imageNoise);

  %--------------------------------------------------------------------------
        % Update Geometry and Number of Frames Based on Displacement Matrix Size
        fprintf('Updating geometry variables based on the size of displacement field matrices.\n')
        [specimen,grid] = func_updateSpecGeom(specimen,grid,disp);

        % Currently the rotations are unused so remove them to save RAM
        disp = rmfield(disp,'rot');
    %--------------------------------------------------------------------------
    % Create the time vector based on the number of frames in the disp struct
    time.numFrames = size(disp.x,3);
    time.vec = 0:time.step:(size(disp.x,3)-1)*time.step;

    %Create arrays of x and y vectors
    X_vec=pos.x;
    Y_vec=pos.y;

 %% Perform displacement corrections
 fprintf('Saving Raw Displacement Fields \n')
            RawDisp=disp;
            fprintf('Correcting Grid Method Displacements along specimen edges\n')
            [disp,DispCorr,grid,ProgramVersions]=func_CorrectGMDisp(disp,...
                DispCorr,grid,pos);

  %% Calculate strain
  % Load the Reference Image and Determine Where the Free Edge is
        fprintf('Obtaining and setting the free edge location.\n')
        [freeEdge,specimen,disp] =... 
        func_getFreeEdge(globalOpts.hardCodeFreeEdge,...
            imagePath,imageFile,specimen,disp);

  % Smooth and Calculate Strain
        fprintf('Calculating strain from the displacement fields.\n')
        [disp,strain,strainRate] = func_smoothCalcStrain(globalOpts,pos,time,...
            grid,disp,smoothingOpts,extrapOpts);

 %% Calculate Acceleration
  % Smooth and Calculate Acceleration
        fprintf('Calculating acceleration from the displacement fields.\n')
        [disp,~,accel] = func_smoothCalcAccel(pos,time,grid,disp,smoothingOpts,...
            extrapOpts,diffOpts);

%% Remove some 3D fields from the structs to save memory.
        if globalOpts.reduceMemory
            disp.x = disp.extrap.x;
            disp.y = disp.extrap.y;
            disp = rmfield(disp,'extrap');
            disp = rmfield(disp,'tSmooth');
        end
%% Save Fields and inputs
fprintf('Saving Kinematic Fields \n')
SaveName=strcat(MainSaveDir,'/',ParentDesig,'_',NoiseCall, ...
    '_KinFields_NoSmooth.mat');
save(SaveName)

%% Record Output Data Structure
OutputStruct.disp=disp;
OutputStruct.RawDisp=RawDisp;
OutputStruct.DispCorr=DispCorr;
OutputStruct.grid=grid;
OutputStruct.specimen=specimen;
OutputStruct.gridMethodOpts=gridMethodOpts;



end