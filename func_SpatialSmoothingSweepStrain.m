function OutputStruct=func_SpatialSmoothingSweepStrain(KernVec, ...
    disp,time,RawStrain,pos,grid, ...
    smoothingOpts,extrapOpts,globalOpts, ...
    SaveDir,Dessig,NoiseCall)

%This function is written to perform a parametric sweep of spatially
    %smoothed displacement fields and to calculate and output an array of
    %the spatially smoothed strain fields

%Author: Andrew Matejunas

%Date Created: 2022-09-22

%Changelog:
    %2022-09-22: Initial version created

%Function Input Variables
    %KernVec- Vector of spatial smoothing kernals
    %disp- structure containing displacement data
    %time- structure containing time information
    %strain- structure containing strain fields
    %pos- structure containing spatial coordinates
    %grid- structure containing grid information
    %smoothgOpts- structure containing general smoothing options
        %NOTE: Existing smoothing kernals must be overwritten
    %extrapOpts- structure containing extrapolation options
    %globalOpts- structure containing global options (unsure if used)
    %SaveDir- Directory in which to save output fields
    %Desig- simulation designation
    %NoiseCall- string indicating whether raw or noisy images are used

%Function Output Variables
    %OutputStruct- Outputs the important variables as a structure
    %SaveFile- Saves the function workspace as a result

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set Base Smoothing Opts
smoothingOpts.FFTempSmooth=false;
smoothingOpts.spatialSmooth=true;

smoothingOpts.spatialFilt='guass';
smoothingOpts.spatialEdgeMode='symmetric';
smoothingOpts.SpaKernVec=KernVec;

%% Initialize Output variables
strain.x=zeros([size(RawStrain.x),length(KernVec)]);
strain.s=zeros([size(RawStrain.x),length(KernVec)]);

%% Record strains with no smoothing
fprintf('Performing smoothing Parametric Sweep \n')
strain.x(:,:,:,1)=RawStrain.x;
strain.s(:,:,:,1)=RawStrain.s;

%% initialize the progress bar
ItTot=length(KernVec);

progNum=1/ItTot;
progMsg=strcat('1/',num2str(ItTot),'smoothing iterations complete');
Progress=waitbar(progNum,progMsg);

%% Perform Smoothing sweep
for k=2:ItTot
    %% set up progress bar
    ItStr=num2str(k);
    progNum=k/ItTot;
    progMsg=strcat(ItStr,'/',num2str(ItTot), ...
        'Smoothing iterations complete');
    %% Update Smoothing Kernal
    smoothingOpts.spatialKernal=KernVec(k)*[1,1];

    %% calculate smoothed strains (ignoring strain rate and disp)
    [~,tempstrain,~]=func_smoothCalcStrain(globalOpts,pos,time,grid,disp, ...
        smoothingOpts,extrapOpts);

    %% Put strains into array
    strain.x(:,:,:,k)=tempstrain.x;
    strain.s(:,:,:,k)=tempstrain.s;
    
   %% Generate Progress bar
   Progress=waitbar(progNum,Progress,progMsg);

end

%% clear unneeded 
clear tempstrain Progress progNum ItStr ItTot progMsg 
%% Save Workspace
fprintf('Saving Smoothing Sweep Data \n')
SaveFile=strcat(SaveDir,'/',Dessig,'_',NoiseCall, ...
    '_SpatialSmoothingSweepStrain.mat');
save(SaveFile)
fprintf('Generating output structure \n')
%% Set Output Structure to output a record of inputs
OutputStruct.smoothingOpts=smoothingOpts;
OutputStruct.extrapOpts=extrapOpts;
OutputStruct.globalOpts=globalOpts;
OutputStruct.disp=disp;
OutputStruct.time=time;
OutputStruct.pos=pos;
OutputStruct.grid=grid;
OutputStruct.File=SaveFile;

%% Output Calculated Strain array
OutputStruct.strain=strain;

end