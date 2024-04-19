function [DispRes,StrainRes,AccelRes,SGRes,StressRes,...
    disp,accel,strain,SG,ConstStress] =...
    func_computeResolutionNoiseFloorV2(NoiseMag,grid,time, ...
    imagePath,imageFile,gridMethodOpts,specimen, ...
    DispCorr,CondOpts,smoothingOpts, ...
    globalOpts,extrapOpts,diffOpts, ...
    material,RefPar)
%This function is written to compute the resolution and noise floor for a
    %grid method sample. Note this code is for synthetic grid images.
    %Experimental images will use a different code

%Author: Andrew Matejunas

%Date created: 2022-12-18

%Version History/change log:
    %V2 2023-01-03: Included calculation of the noise floor for
        %constitutive model stresses assuming the constituitve model
        %parameters are known
    %2023-01-04: Changed code to incorporate different smoothing parameters
    %2023-01-08: Added calculation of the noise floor in the average strain
        %along a vertical slice some distance X from the free surface of an
        %IBII specimen. This is to better allow thresholding the cost
        %function to account for the noise floor. 

%Function Input Parameters
    %NoiseMag- Magnitude (% of dynamic range)
    %grid- structure containing grid method parameters
    %time-structure containing time information
    %imagePath- folder containing the undeformed grid images
    %imageFile- filename for the first image in the sequence
    %gridMethodOpts- struct containing grid method processing options
    %specimen- struct containing specimen information
    %DispCorr- struct containing kinematic field correction options
    %CondOpts- struct containing conditioning data
    %material- structure containing material information
    %RefPar- Structure containing reference consitutive parameters
    %smooothingOpts- Structure containing the the spatial and temporal
        %smoothing parameters to use when calculating strain and 
        %acceleration




%Function Output Arguments
   %DispRes- Struct containing Displacentent Resolution & Noise floor in 
        %avg- mean
        %NF- Resolution/Noise Floor
   %StrainRes- Struct containing strain Resolution & Noise floor in 
        %av- mean
        %NF-  Resolution/Noise Floor
   %AccelRes- Struct containing acceleration Resolution & Noise floor in 
        %avg- mean
        %NF-  Resolution/Noise Floor
   %SGRes- Struct containing Resolution & Noise floor in Stress guage 
        %stresses
   %Stress Res- structure constaining the resolution in noise floor in
        %stresses calculated using the reference constitutive parameters
   %disp- struct containing full field displacements
   %accel- struct containing full field accelerations
   %strain- struct containing full-field in-plane strains
   %SG- structure containing stress gauge
   %ConstStress- stresses calculated with the constitutive model and the
        %reference parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define Noise Parameters
imageNoise.addNoise=1;
imageNoise.pcNoise=NoiseMag;
imageNoise.bits=16;
imageNoise.convToUInt16=1;

%% Remove some options from original IBII codes

FEValidMode=0;
calcKinFieldsFromDisp=1;



%% %Preallocate identified constitutive parameter
% %     fprintf(strcat('Processing iteration',itNum,' \n'))
%     fprintf('------------------------------------------------------------- \n')
    %% IMAGE PROCESSING: Use the Grid Method to extract displacement fields
    if ~FEValidMode
        % Process the raw tiff images using the grid method code developed by
        % Grediac et al.
        fprintf('\n--------------------------------------------------------------\n')
        fprintf('GRID METHOD PROCESSING\n')
        fprintf('--------------------------------------------------------------\n')

        %--------------------------------------------------------------------------
        % GRID IMAGE PROCESSING



              
%             fprintf('Processing images using the grid method toolbox.\n')
            % Process the image squence with the grid method toolbox
            [grid,pos,disp] = func_gridMethodImageProcessing_AJM(imagePath,...
                imageFile,...
                grid,gridMethodOpts,imageNoise);    
             fprintf('Grid Method Processing Complete.\n')
       

        %--------------------------------------------------------------------------
        % Update Geometry and Number of Frames Based on Displacement Matrix Size
%         fprintf('Updating geometry variables based on the size of displacement field matrices.\n')
        [specimen,grid] = func_updateSpecGeom(specimen,grid,disp);

        % Currently the rotations are unused so remove them to save RAM
        disp = rmfield(disp,'rot');
    end

    %--------------------------------------------------------------------------
    % Create the time vector based on the number of frames in the disp struct
    time.numFrames = size(disp.x,3);
    time.vec = 0:time.step:(size(disp.x,3)-1)*time.step;

    %Create arrays of x and y vectors
    X_vec=pos.x;
    %Y_vec=pos.y;

    %% POST-PROCESSING: Smoothing and Kinematic Field Derivation
    % Smooth the displacement data and then calculate acceleration and strains
    % Extrapolate the data to account for the missing pitch on the edges
    fprintf('\n--------------------------------------------------------------\n')
    fprintf('POST PROCESSING: Smoothing and Kinematic Field Calculation\n')
    fprintf('--------------------------------------------------------------\n')
    %% Perfom corrections if needed
%     
%             fprintf('Saving Raw Displacement Fields \n')
%             RawDisp=disp;
            fprintf('Correcting Grid Method Displacements along specimen edges')
            [disp,DispCorr,grid,~]=func_CorrectGMDisp(disp,...
                DispCorr,grid,pos);
            
       


    % Calculate the kinematic fields from displacement fields using
    % displacements from images or displacements from FE data
    if calcKinFieldsFromDisp
        %--------------------------------------------------------------------------
        % Load the Reference Image and Determine Where the Free Edge is
%         fprintf('Obtaining and setting the free edge location.\n')
        [~,~,disp] = func_getFreeEdge(globalOpts.hardCodeFreeEdge,...
            imagePath,imageFile,specimen,disp);

        %--------------------------------------------------------------------------
        % Smooth and Calculate Strain
        fprintf('Calculating strain from the displacement fields.\n')
        [disp,strain,~] = func_smoothCalcStrain(globalOpts,pos,time,...
            grid,disp,smoothingOpts,extrapOpts);
        %Post Correct Strain
       strain=func_PostCorrectStrain(grid,strain,pos,time,DispCorr,...
                    smoothingOpts);
        %--------------------------------------------------------------------------
        % Smooth and Calculate Acceleration
        fprintf('Calculating acceleration from the displacement fields.\n')
        [disp,~,accel] = func_smoothCalcAccel(pos,time,grid,disp,smoothingOpts,...
            extrapOpts,diffOpts);

        % Remove some 3D fields from the structs to save memory.
        if globalOpts.reduceMemory
            disp.x = disp.extrap.x;
            disp.y = disp.extrap.y;
            disp = rmfield(disp,'extrap');
            disp = rmfield(disp,'tSmooth');
        end
%     else
%         fprintf('Using kinematic fields directly from FE data.\n')
%         % For pure FE data the kinematic fields can be taken directly from the
%         % FE calculation
%         disp.xAvg = func_avgFFVarOverWidth(disp.x);
%         accel.xAvg = func_avgFFVarOverWidth(accel.x);
%         strain.xAvg = func_avgFFVarOverWidth(strain.x);
%         strain.yAvg = func_avgFFVarOverWidth(strain.y);
%         [strainRate,~] = func_calculateStrainRate(strain,time,smoothingOpts,true);
%         strainRate.xAvg = func_avgFFVarOverWidth(strainRate.x);
%         strainRate = func_calcMaxStrainRate(smoothingOpts,grid,strainRate);
%         disp.extrap.x = disp.x;
%         disp.extrap.y = disp.y;
    end
    %% Calculate Stress Gage Stresses

    fprintf('Calculating Stress Gage Stresses \n')
    SG=func_Full_SG(accel,X_vec,time,material.rho);

  %% Condition the data (censoring/downsampling/smoothing)
    %% Define Censorship Parameter. (Standard is 4 grid pitches from
        %impact and free edges)
        fprintf('Conditioning IBII Fields \n')
        [SG,accel,strain,~]=func_conditionIBIIData(SG,accel,...
            strain,X_vec,time,CondOpts);

%% Calculate Constitutive Model stresses
ConstStress=func_ViscoConstitutiveV6(strain.x,strain.y,strain.s,time.vec, ...
    RefPar,0,0,0);
%% Calculate average strains
AvXstrain=squeeze(mean(strain.x));
AvYstrain=squeeze(mean(strain.y));
AvSstrain=squeeze(mean(strain.s));

%% Calculate resolution and Noise floor
    %Resolution is the same as noise floor and is the standard deviation of
    %the full field data
    
dispStdx= std(disp.x,0,"all");
dispAvgx=mean(disp.x,"all");
DispRes.x.NF=dispStdx;
DispRes.x.avg=dispAvgx;

dispStdy = std(disp.y,0,"all");
dispAvgy=mean(disp.y,"all");
DispRes.y.NF=dispStdy;
DispRes.y.avg=dispAvgy;

strainStdx=std(strain.x,0,"all");
strainAvgx=mean(strain.x,'all');
StrainRes.x.NF=strainStdx;
StrainRes.x.avg=strainAvgx;
StrainRes.xAv.NF=std(AvXstrain,0,"all");
StrainRes.xAv.avg=mean(AvXstrain,"all");


strainStdy=std(strain.y,0,"all");
strainAvgy=mean(strain.y,'all');
StrainRes.y.NF=strainStdy;
StrainRes.y.avg=strainAvgy;
StrainRes.yAv.NF=std(AvYstrain,0,"all");
StrainRes.yAv.avg=mean(AvYstrain,"all");


strainStds=std(strain.s,0,"all");
strainAvgs=mean(strain.s,'all');
StrainRes.s.NF=strainStds;
StrainRes.s.avg=strainAvgs;
StrainRes.sAv.NF=std(AvSstrain,0,"all");
StrainRes.sAv.avg=mean(AvSstrain,'all');

accelStdx=std(accel.x,0,"all");
accelAvgx=mean(accel.x,'all');
AccelRes.x.NF=accelStdx;
AccelRes.x.avg=accelAvgx;

accelStdy=std(accel.y,0,"all");
accelAvgy=mean(accel.y,'all');
AccelRes.y.NF=accelStdy;
AccelRes.y.avg=accelAvgy;

SGstdx=std(SG.x,0,"all");
SGavgx=mean(SG.x,'all');
SGRes.x.NF=SGstdx;
SGRes.x.avg=SGavgx;

SGstds=std(SG.s,0,"all");
SGavgs=mean(SG.s,'all');
SGRes.s.NF=SGstds;
SGRes.s.avg=SGavgs;

StressSDx=std(ConstStress.xx,0,'all');
Stressavgx=mean(ConstStress.xx,'all');
StressRes.x.NF=StressSDx;
StressRes.x.avg=Stressavgx;

StressSDy=std(ConstStress.yy,0,'all');
Stressavgy=mean(ConstStress.yy,'all');
StressRes.y.NF=StressSDy;
StressRes.y.avg=Stressavgy;

StressSDs=std(ConstStress.xy,0,'all');
Stressavgs=mean(ConstStress.xy,'all');
StressRes.s.NF=StressSDs;
StressRes.s.avg=Stressavgs;

end
