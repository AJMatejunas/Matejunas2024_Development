function OutputStruct=func_SmoothSaveKinFields(globalOpts, ...
            extrapOpts,smoothingOpts,diffOpts, ...
            pos,time,grid,disp, ...
            MainSaveDir,ParentDesig,NoiseCall,SmoothLevel, ...
            DispCorr,specimen,imageNoise) 

%This function is intended to smooth displacement fields, and calculate
    %strain and accelerations fields. The code also saves the output
    %kinematic fields into the working directory for future use.

    %Author: Andrew Matejunas
    %Date Created: 2022-09-19
    
    %Changelog/Version History:
        %2022-09-19- Initial function written
        %2022-09-21- Added capability to calculate average X and Shear
            %strains

% Function Input Parameters:
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
       %globalOpts- structure containing more options
    %extrapOpts- structure containing edge extrapolation options
    %diffOpts- structure containing temporal differentiation options
    %SmoothLevel- String identifying amount of smoothing performed
    %smoothingOpts- Structure describing the smoothing operations to be
        %performed
    %disp- Unsmoothed displacement fields
    %DispCorr- Record of displacement correction parameters
    %specimen- Record of specimen details
    %imageNoise- record of 

% Function outputs
    %OuputStruct- Structure containing output kinematic fields 
    %SaveName- Saved file in the save directory


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Save Displacements Before smoothing
UnsmoothDisp=disp;
X_vec=pos.x;
Y_vec=pos.y;

%% Smooth and Calculate Strain
[disp,strain,strainRate] = func_smoothCalcStrain(globalOpts,pos,time,...
            grid,disp,smoothingOpts,extrapOpts);

%% Smooth and Calculate Accelerations
[disp,~,accel] = func_smoothCalcAccel(pos,time,grid,disp,smoothingOpts,...
            extrapOpts,diffOpts);

%% Calculate average strains
avgX_strain=squeeze(mean(strain.x,1));
avgXY_strain=squeeze(mean(strain.s,1));


%% Save Kinematic fields
SaveFile=strcat(MainSaveDir,'/',ParentDesig,'_',NoiseCall,'_', ...
    SmoothLevel,'Smooth.mat');
save(SaveFile)
clear MainSaveDir

%% Output Kinematic Fields
OutputStruct.UnsmoothDisp=UnsmoothDisp;
OutputStruct.disp=disp;
OutputStruct.strain=strain;
OutputStruct.strainRate=strainRate;
OutputStruct.accel=accel;
OutputStruct.pos=pos;
OutputStruct.time=time;
OutputStruct.X_vec=X_vec;
OutputStruct.Y_vec=Y_vec;
OutputStruct.avgX_strain=avgX_strain;
OutputStruct.avgXY_strain=avgXY_strain;


%% Ensure Input Parameters are in the output
OutputStruct.globalOpts=globalOpts;
OutputStruct.extrapOpts=extrapOpts;
OutputStruct.smoothingOPpts=smoothingOpts;
OutputStruct.diffOpts=diffOpts;
OutputStruct.grid=grid;
OutputStruct.ParentDesig=ParentDesig;
OutputStruct.File=SaveFile;
OutputStruct.imageNoise=imageNoise;
OutputStruct.DispCorr=DispCorr;
OutputStruct.specimen=specimen;
OutputStruct.NoiseCall=NoiseCall;
OutputStruct.SmoothLevel=SmoothLevel;




end