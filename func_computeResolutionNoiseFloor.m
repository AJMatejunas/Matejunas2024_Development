function [DispRes,StrainRes,AccelRes,SGRes,disp,accel,strain,SG] =...
    func_computeResolutionNoiseFloor(NoiseMag,grid,time,imagePath,...
    imageFile,gridMethodOpts,specimen,DispCorr,CondOpts,globalOpts, ...
    extrapOpts,diffOpts,material)
%This function is written to compute the resolution and noise floor for a
    %grid method sample. Note this code is for synthetic grid images.
    %Experimental images will use a different code

%Author: Andrew Matejunas

%Date created: 2022-12-18

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
    %SGRes- Struct containing Displacentent Resolution & Noise floor in 
   %disp- struct containing full field displacements
   %accel- struct containing full field accelerations
   %strain- struct containing full-field in-plane strains
   %SG- structure containing stress gauge




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create smoothing opts data structure without smoothing
smoothingOpts.spatialSmooth=0;
smoothingOpts.WATempSmooth=0;
smoothingOpts.FFTempSmooth=0;
smoothingOpts.spatialFilt='gauss';
smoothingOpts.spatialKernal=[3,3];
smoothingOpts.spatialEdgeMode='symmetric';
smoothingOpts.FFTemporalFilt='sgolay';
smoothingOpts.FFTemporalPad=0;
smoothingOpts.FFTemporalPadFrames=3;
smoothingOpts.FFTemporalPadMethod='replicate';
smoothingOpts.WATemporalAvgFirst=0;
smoothingOpts.WATemporalFilt='sgolay';
smoothingOpts.WATemporalKernal=[11,3];
smoothingOpts.WATemporalPad=0;
smoothingOpts.WATemporalPadFrames=3;
smoothingOpts.WATemporalPadMethod='replicate';

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
  
%% Calculate Displacement resolution and Noise floor
    %Resolution is
    
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

strainStdy=std(strain.y,0,"all");
strainAvgy=mean(strain.y,'all');
StrainRes.y.NF=strainStdy;
StrainRes.y.avg=strainAvgy;

strainStds=std(strain.s,0,"all");
strainAvgs=mean(strain.s,'all');
StrainRes.s.NF=strainStds;
StrainRes.s.avg=strainAvgs;

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
SGavgs=mean(SG.s,'all')
SGRes.s.NF=SGstds;
SGRes.s.avg=SGavgs;

end
