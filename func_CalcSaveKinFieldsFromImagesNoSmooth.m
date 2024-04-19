function OutputStruct =...
    func_CalcSaveKinFieldsFromImagesNoSmooth(MainSaveDir,...
    imageFile,imagePath,ParentDesig,NoiseCall, ...
        grid,specimen,time, ...
        gridMethodOpts,imageNoise,DispCorr,globalOpts,extrapOpts,...
       smoothingOpts,diffOpts)
%This function is written to calculate and save the displacement, strain,
    %and acceleration kinematic fields from image data and processing
    %parameters.

    %Author: Andrew Matejunas
    %Date Completed: 2022-09-19

    %Change Log/Version History:
        %2022-09-19: Code created. First iteration
        %2022-09-21: Fixed smoothingOpts bug

% Function Input Parameters
    %MainSaveDir- Directory to save kinematic fields in
    %imageFile- Filename of the first grid image in the sequence
    %imagePath- Directory containing grid images
    %ParentDesig- Designation of the FE file used to generate grid images
        %or the experiment from which the images were produced
    %NoiseCall- String identifying whether or not to add noise to the
        %displacement fields
    %grid- Structure containing grid method parameters
    %specimen- Structure containg specimen dimensions
    %time- structure containing time information
    %gridMethodOpts- structure containing grid method processing options
    %imageNoise- Structure describing how to add noise to the raw images
    %DispCorr- structure containing the information on how to correct for
        %corrupted displacement data on the grid edges
    %globalOpts- structure containing more options
    %extrapOpts- structure containing edge extrapolation options
    %diffOpts- structure containing temporal differentiation options

% Function Outputs
    %OutputStruct- Structure containing a record of input parameters and
        %the generated displacement, strain, and acceleration fields
    %SaveName- Kinematic fields and inputs saved to a directory on the
        %computer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set Smoothing options
%Hardcode No smoothing
smoothingOpts.spatialSmooth=0;
smoothingOpts.WATempSmooth=0;
smoothingOpts.FFTempSmooth=0;
smoothingOpts.spatialKernal=[0,0];
smoothingOpts.FFTemporalKernal=[0,3];


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
%% Calculate average strains
avgX_strain=squeeze(mean(strain.x,1));
avgXY_strain=squeeze(mean(strain.s,1));

%% Save Fields and inputs
fprintf('Saving Kinematic Fields \n')
SaveName=strcat(MainSaveDir,'/',ParentDesig,'_',NoiseCall, ...
    '_KinFields_NoSmooth.mat');
clear MainSaveDir
save(SaveName)


%% Record Output Data Structure
OutputStruct.disp=disp;
OutputStruct.time=time;
OutputStruct.RawDisp=RawDisp;
OutputStruct.DispCorr=DispCorr;
OutputStruct.grid=grid;
OutputStruct.specimen=specimen;
OutputStruct.gridMethodOpts=gridMethodOpts;
OutputStruct.imageNoise=imageNoise;
OutputStruct.globalOpts=globalOpts;
OutputStruct.pos=pos;
OutputStruct.X_vec=X_vec;
OutputStruct.Y_vec=Y_vec;
OutputStruct.NoiseCall=NoiseCall;
OutputStruct.ProgramVersions=ProgramVersions;
OutputStruct.freeEdge=freeEdge;
OutputStruct.strain=strain;
OutputStruct.strainRate=strainRate;
OutputStruct.accel=accel;
OutputStruct.File=SaveName;
OutputStruct.extrapOpts=extrapOpts;
OutputStruct.diffOpts=diffOpts;
OutputStruct.avgX_strain=avgX_strain;
OutputStruct.avgXY_strain=avgXY_strain;


end