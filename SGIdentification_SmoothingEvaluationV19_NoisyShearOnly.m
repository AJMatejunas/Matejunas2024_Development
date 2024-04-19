%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NOTE: THIS SCRIPT IS ONLY FOR USE WITH GRID IMAGES. USE OTHER CODES FOR7
%FINITE ELEMENT VALIDATION.


%This code is written to evaluate the effects of smoothing and noise on the
%identification of viscoelastic constitutitve parameters using the image
%based inertial impact test. 

%Authors: Andrew Matejunas, 
% THIS CODE IS HEAVILY ADAPTED FROM main_IBIIProcessing_v1_0r 
%Date completed: 2022-06-24

%Version History/changelog
 %V2- added 3 initial guesses to check convergence
 %V3- added correction to shear and Y strain
   %2022-11-22: Added plotting in a separate function in attempt to prevent
    %GPU crashes
 %V4- 2022-12-8 Perfomance increase pass to improve computational time
      %Narrowed upper and lower bounds on parameter identification
      %increased step tolerance for fmincon to 1e-3
      %increased downsampling in X to 6 pixels so that downsampling is
        %comparible
 %V5- 2022-12-9: Made spacial smoothing Kernel asymmetric
      %2022-12-13: Reverted to symmetric Smoothing Kernels
      %2023/01/12: Converted the cost function to only consider and
        %identify shear
      %2023/01/18: Changed script file to SGIdentification_
                   %SmoothingEvaluationV5_SeparateKG_SameSmooth. Added in
                   %the capability to identify bulk modulus after first
                   %identifying shear modulus and the associated time
                   %constant
 %V6- 2023/02/09: Changed the strain correction algorithm to 
                  %func_PostCorrectStrainV2 to remove the erroneous
                  %strainyy=0 boundary condtion at the free surface and
                  %schanges the corection of normal strain in the y
                  %direction to a simple linear extrapolation
      %2023/02/13: Changed the correction interval on strain and
                     %displacement to be equal at 1.4P
                   %Narrowed the bounds on the initial guess and minimized
                   %values to +\- 25%
 %V7- 2023/02/16:  Changed displacement correction to a quadradic fit
                        %(func_CorrectGMDispV3) and a further change to the
                        %strain correction (func_PostCorrectStrainV3)
 %V8- 2023/02/16: Added a linear correction to the normal strains in the x
                        %direction and changed the strain correction
                        %function to func_PostCorrectStrainV3
 %V9- 2023/02/22: Added the ability to censor frames from the end of the
                    %Data set (func_conditionIBIIDataV2) and also added
                    %corrections to the acceleration fields on all 4 edges
                    %(func_PostCorrectFields)
%V10- 2023/02/23: Removed correction on the x strains and accelerations on
                    %the top and bottom edges
                    %(func_PostCorrectFieldsNoTBX) Also changed the free
                    %end censorship back to 5 grid pitches
%V11- 2023/02/27:  Recoupled K and G to bring back the original cost
                     %function with the improved edge correction algorithm
                     %and renamed program to 
                     %SGIdentification_SmoothingEvaluationV11_KG_Coupled
                     %using 
                   %Changed edge correction back to func_PostCorrectFields
                   %Added outputs of phi_K and phi_G
%V12- 2023/03/01: Changed version of minimization algorithm  
                  %func_PronySGvfmV9 allowing for weighting of cost
                  %function
%V12noShear- weighted phi_G at 0% weighting
%V13- 2023/03/02- changed displacement correction algorithim to
                    %func_CropAndExtrapolateDisp_V4
                  %Changed strain calculation, smoothing, and correction
                    %algorithm to func_smoothCalcStrain_v4
                  %Changed acceleration calculation, smoothing, and
                    %correction algorithm to func_smoothCalcAccel_v4.
                  %These changes also involved a large overhaul of the
                    %processing parameters data stuctures
                  %changed minimization algorithm to func_PronySGvfmV10
                    %allowing for custom weighting of the cost function
%V14- 2023/03/05: Changed code to obtain CondOpts, SolveOpts, and minOpts
                        %structures from the processing parameters data 
                        %file. Removed options to specify them from the 
                        %code 
%V15_ShearOnly- 2023/03/13: Remove identification of bulk modulus
                           %Updated strain extrapolation pramater to
                            %extrapolate linearly over 1/2 grid pitch +1 px
                            %or 7px whichever is greater
                           %Fixed quadratic extrapolation of accelerations
                            %to 7px for accelerations
%V15_BulkOnly-   2023/03/16: Changed from identification of shear modulus
                                %to identification of bulk modulus for both
                                %exact values of shear modulus and the
                                %shear modulus identified from the pure
                                %shear cost function
%V15_NoisyShearOnly- 2023/03/16: Changed the way noise is handled to keep up
                                %with the noise free smoothing sweep
                                %changes and added noise
%V16_NoisyShearOnly- 2023/03/16: Now solve for noise iterations within a
                                    %parfor loop. Significant processing 
                                    %time improvements
%V17_NoisyShearOnly- 2023/03/17: Turns off temporal padding. Instead start
                                    %cropping 1/2 temporal smoothing kernal 
                                    %from the start and end of the image set 
%V19_NoisyShearOnly- 2023/03/27: Made spatial cropping and extrapolation
                                       %intervals the same
                     %2023/08/17: Removed custom color map option. For
                     %perceptually uniform color maps, use crameri.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize (clear memory and figures)
close all
clear variables
clc



%% Remove some options from original IBII codes
hardCodePath=0;
FEValidMode=0;
calcKinFieldsFromDisp=1;
saveGridData=0;

%% Create reasonable upper bounds
%bounds on the initial guess are +/- 25%
%% Create reasonable upper bounds
%bounds on the initial guess are +/- 25%
ub=[1.875e9;
    1.0962E9;... G1
    12.5E-6]; %tau_1
ubG=ub(2:3);

%% Create reasonable lower bounds
lb=[1.125e9;
    657.7E6;... G1
    7.5E-6]; %tau_1
lbG=lb(2:3);

%% Create Matrix of initial guesses
    %Eint=2.5E9;
    %Choose a random 3x3 matrix of initial guress multipliers. Multipliers
        %will be uniformly distributed random numbers in the interval [0,1]
 
    intMult=rand(3,3);
    intMultK=rand(1,3);
  
    %Create 3X3 matrix of initial guess parameters
     intMatrix=lb+(ub-lb).*intMult;
     intMatrix=intMatrix(2:3,:);

%     intGuess=[Eint/(3*(1-2*nu));... K1
%               Eint/(2*(1+nu));.... G1
%               10E-6];... tau1
   

%% Add reference parameters file 
[RefPar.file,RefPar.path]=uigetfile('*.mat',...
    'Choose Mat file containing reference constitutive pamrameters');
RefPar=load(strcat(RefPar.path,'/',RefPar.file),'MatProps');
%% Create structute of known parameters
    Einf=2.98E9;
    nu=0.26;

    exactProps.Kinf=Einf/(3*(1-2*nu));
    exactProps.Ginf=Einf/(2*(1+nu));
    exactProps.nu=0;
    exactProps.Ki=RefPar.MatProps.Ki;
   
 
%% Choose Directory to save results
savePath=uigetdir({},'Choose Folder to Save Results');

%% INITIALISE: Add path for processing functions
%(this section taken from main_IBIIProcessing_v1_0r
fprintf('Adding pathes for processing functions.\n')
% Add the path for the grid method code and other useful functions:
funcPath = [pwd,'\Functions\'];

% If the default path is not found we should find it
if exist(funcPath,'file') ~= 7
    hWarn = warndlg('Folder for processing functions not found.',...
        'Function folder not found');
    waitfor(hWarn);
    funcPath = uigetdir(pwd,'Locate Processing Function Folder');
end
addpath(funcPath);
addpath([funcPath,'GridMethodToolbox\']);

%% LOAD RAW DATA FILES: Raw .tiff images
    
    fprintf('Loading reference image from the selected test data folder.\n')
    if ~hardCodePath
        [imageFile,imagePath] = uigetfile({'*.*','All Files'},'Select the first image in the sequence');
    else   
        imageFile = 'DefGridImage_001.tiff';
        imagePath = 'E:\Matlab_WorkingDirectory\1_IBIITest_Data\TestData_ID_CFIPQIso_HalfPulse\'; 
    end

%% Load Processing Parameter Data Structures
fprintf('Loading processing parameters file.\n')
initPath = imagePath;
initFile = 'processingParameters.mat';
if exist([initPath,initFile],'file') ~= 2
    hWarn = warndlg('Processing parameter file does not exist.','Processing parameters not found');
    waitfor(hWarn);
    [initFile,initPath,~] = uigetfile('*.mat','Locate processing parameter file');
end
% Load the processing parameters from file


load([initPath,initFile])

% Store the FE valid mode parametervccvv.MatProps.Gi*100;
globalOpts.FEValidMode = FEValidMode;
globalOpts.calcKinFieldsFromDisp = calcKinFieldsFromDisp;

%% Define Test Dseignation
TestDeg=char(cell2mat(inputdlg('Input Test Designation')));

 
%% Define Smoothing Opts
TempKernVec=[0,5:2:31];
%TempKernVec=0;
TempKernNum=length(TempKernVec);  
%SpaKernVec=[0,5:2:21];
SpaKernVec=43:2:51;
%SpaKernVec=0;

%record Number of Kernels
SpaKernNum=length(SpaKernVec);

%% Define whether Noise should be added
imageNoise.addNoise=true;
imageNoise.numCopies=30;
NumNoise=imageNoise.numCopies;
 

%% Generate Progress bar
itCount=0;
itTotal=TempKernNum*SpaKernNum;

progmsg=strcat('0/',num2str(itTotal),...
    ' Smoothing Iterations Complete');
Progress=waitbar(0,progmsg);

%% set the extrapolation options that do not depend on 
extrapOpts.accel.tempPadOn=false;       


[strainCorr,strainSmooth,dispCorr,accelCorr]=deal(cell(SpaKernNum,1));

for m=1:length(SpaKernVec)
            if SpaKernVec(m)==0
            smoothOpts.strain.spatialSmooth=false;
            % Set default displacement cropping and extrapolation
            %parameters at 7 pixels
            extrapOpts.disp.extrapPx1st=7;
            extrapOpts.disp.extrapPx2nd=7;
            extrapOpts.disp.cropPx1st=7;
            extrapOpts.disp.cropPx2nd=7;
        else
            smoothOpts.strain.spatialSmooth=true;
            smoothOpts.strain.spatialKernelSize=SpaKernVec(m)*[1,1];
            smoothOpts.strain.spatialKernelStd=[10,10];

            % Scale displacement cropping and extrapolation parameters
                %with smoothing Kernel size
                if SpaKernVec(m)<=13
                   extrapOpts.disp.cropPx1st=7;
                   extrapOpts.disp.cropPx2nd=7;
                   extrapOpts.disp.extrapPx1st=extrapOpts.disp.cropPx1st;
                   extrapOpts.disp.extrapPx2nd=extrapOpts.disp.cropPx2nd;
                else

                   extrapOpts.disp.cropPx1st=round(SpaKernVec(m)/2)+1;
                   extrapOpts.disp.cropPx2nd=round(SpaKernVec(m)/2)+1;
                   extrapOpts.disp.extrapPx1st=extrapOpts.disp.cropPx1st;
                   extrapOpts.disp.extrapPx2nd=extrapOpts.disp.cropPx2nd;

                end
        
         % Adjust sample conditioning to remove at least one extrapolation
            %kernel from the cost function evaluation

        if extrapOpts.strain.extrapPx1st<=20
            CondOpts.FreeCens=20;
        else
            CondOpts.FreeCens=extrapOpts.strain.extrapPx1st;
        end

            end

        extrapOpts.strain.extrapPx1st=extrapOpts.disp.extrapPx1st;
        extrapOpts.strain.extrapPx2nd=extrapOpts.disp.extrapPx2nd;
        extrapOpts.strain.enforceGlobBCsPx=...
            extrapOpts.strain.extrapPx1st;
        extrapOpts.accel.extrapPx1st=extrapOpts.disp.extrapPx1st;
        etrapOpts.accel.extrapPx2nd=extrapOpts.disp.extrapPx2nd;

        accelCorr{m}=extrapOpts.accel;
        strainSmooth{m}=smoothOpts.strain;
        strainCorr{m}=extrapOpts.strain;
        dispCorr{m}=extrapOpts.disp;
end

%% Initialize contitutive Parameter 
ConstitutiveParam=zeros(TempKernNum,SpaKernNum,NumNoise,2);
Rho=material.rho;
timefull=time;
globalEdge=globalOpts.hardCodeFreeEdge;
for k=1:TempKernNum        
    
    if TempKernVec(k)==0
        smoothOpts.accel.temporalSmooth=false;

    else
        smoothOpts.accel.temporalSmooth=true;
        smoothOpts.accel.temporalKernelSize=TempKernVec(k)*[1,1];
        smoothOpts.accel.temporalKernelOrder=[3,3];
    end   

    if smoothOpts.accel.temporalSmooth(1)==true
           CondOpts.cutStartFrames=...
                round(smoothOpts.accel.temporalKernelSize(1)/2);
           if CondOpts.cutStartFrames<=4
               CondOpts.cutEndFrames=4;
           else
               CondOpts.CutEndFrames=CondOpts.cutStartFrames;
           end
    end

    accelSmooth=smoothOpts.accel;
    for m=1:SpaKernNum
         itCount=itCount+1;
         itNum=num2str(itCount);

        CorrDisp=cell2mat(dispCorr(m));
        SmoothStrain=cell2mat(strainSmooth(m));
        CorrStrain=cell2mat(strainCorr(m));
        CorrAccel=cell2mat(accelCorr(m));
             
        parfor n=1:NumNoise
            time=timefull;
%             fprintf('------------------------------------------------------------- \n')
%             fprintf(strcat('Processing iteration',itNum,' \n'))
%             fprintf('------------------------------------------------------------- \n')
            % IMAGE PROCESSING: Use the Grid Method to extract displacement fields
            
                % Process the raw tiff images using the grid method code developed by
                % Grediac et al.
%                 fprintf('\n--------------------------------------------------------------\n')
%                 fprintf('GRID METHOD PROCESSING\n')
%                 fprintf('--------------------------------------------------------------\n')
% 
%                 %--------------------------------------------------------------------------
%                 % GRID IMAGE PROCESSING
                 
                   % fprintf('Processing images using the grid method toolbox.\n')
                    % Process the image squence with the grid method toolbox
                    [gridN,pos,disp] = func_gridMethodImageProcessing_AJM(imagePath,...
                        imageFile,...
                       grid,gridMethodOpts,imageNoise);                                 
                                      
                %--------------------------------------------------------------------------
                % Update Geometry and Number of Frames Based on Displacement Matrix Size
                %fprintf('Updating geometry variables based on the size of displacement field matrices.\n')
                [specimenN,gridN] = func_updateSpecGeom(specimen,gridN,disp);

                % Currently the rotations are unused so remove them to save RAM
                %disp = rmfield(disp,'rot');
           

            %--------------------------------------------------------------------------
            % Create the time vector based on the number of frames in the disp struct
           % time.numFrames = size(disp.x,3);
            %time.vec = 0:time.step:(size(disp.x,3)-1)*time.step;

            %Create arrays of x and y vectors
            X_vec=pos.x;
            Y_vec=pos.y;
            pos.lengthX = pos.x(end)+pos.xStep/2;
            pos.lengthY = pos.y(end)+pos.yStep/2;
            pos.xGridF = padarray(pos.xGrid,[0,0,time.numFrames-1], ...
                'replicate','post');
            pos.yGridF = padarray(pos.yGrid,[0,0,time.numFrames-1], ...
                'replicate','post');
            pos.x0F = squeeze(padarray(pos.x,[0,0,time.numFrames-1], ...
                'replicate','post'));

            % POST-PROCESSING: Smoothing and Kinematic Field Derivation
            % Smooth the displacement data and then calculate acceleration and strains
            % Extrapolate the data to account for the missing pitch on the edges
%             fprintf('\n--------------------------------------------------------------\n')
%             fprintf('POST PROCESSING: Smoothing and Kinematic Field Calculation\n')
%             fprintf('--------------------------------------------------------------\n')
            % Crop and Extrapolate displacement fields
            [disp.x,disp.y,disp.rX,disp.rY]=...
                func_cropAndExtrapDispFields_v4(pos,disp.x,disp.y, ...
                CorrDisp,false);
            

            % Calculate the kinematic fields from displacement fields using
            % displacements from images or displacements from FE data
           
                %--------------------------------------------------------------------------
                % Load the Reference Image and Determine Where the Free Edge is
                %fprintf('Obtaining and setting the free edge location.\n')
                [freeEdge,specimenN,disp] =...
                    func_getFreeEdge(globalEdge,...
                    imagePath,imageFile,specimenN,disp);

                %--------------------------------------------------------------------------
                % Smooth and Calculate Strain
               % fprintf('Calculating strain from the displacement fields.\n')

                [strain,~]=func_smoothCalcStrain_v4(pos,time,disp,...
                    SmoothStrain,CorrStrain,false);

                %--------------------------------------------------------------------------
                % Smooth and Calculate Acceleration
                %fprintf('Calculating acceleration from the displacement fields.\n')
                [accel,~,~] = func_smoothCalcAccel_v4(pos,time,disp, ...
                    accelSmooth,...
                    CorrAccel,diffOpts,false);

                % Remove some 3D fields from the structs to save memory.
               
            % Calculate Stress Gage Stresses

            %fprintf('Calculating Stress Gage Stresses \n')
            SG=func_Full_SG(accel,X_vec,time,Rho);

            % Condition the data (censoring/downsampling/smoothing)
            %fprintf('Conditioning IBII Fields \n')
            [SG,accel,strain,X_vec,time]=...
            func_conditionIBIIDataV2(SG,accel,...
                strain,X_vec,time,CondOpts);

            % Run the minimization
            [TempG,TempTau,TempPhi]=deal(zeros(3,1));
            for ig=1:3
              intGuess=squeeze(intMatrix(:,ig)); %#ok<PFBNS> 
                      
                [TempParam,TempPhi(ig)]=func_PronyShearSGvfm(strain.s, ...
                    time.vec,SG.s,exactProps,...
                    intGuess,ubG,lbG,SolveOpts,minOpts,CondOpts);
                    TempG(ig)=TempParam(1);
                    TempTau(ig)=TempParam(2);

               
            end
           
                minPhi=min(TempPhi);
                G=TempG(TempPhi==minPhi);
                tau=TempTau(TempPhi==minPhi);
                phi(k,m,n,:)=TempPhi';
                ConstitutiveParam(k,m,n,:)=[G(1);tau(1)];
            
        end
   
                progmsg=strcat(itNum,'/',num2str(itTotal),...
                    ' Smoothing Iterations complete');
                Progress=waitbar(itCount/itTotal,Progress,progmsg);
   %free up memory
   clear disp strain accel SG gridN
    end
end
close(Progress)

%% Separate out the identified parameters
Ident.G=squeeze(ConstitutiveParam(:,:,:,1));
Ident.tau=squeeze(ConstitutiveParam(:,:,:,2));

%% Calculate Identification Errors
Errors.G=(Ident.G...
    -RefPar.MatProps.Gi)/RefPar.MatProps.Gi*100;

Errors.tau=(Ident.tau...
    -RefPar.MatProps.tau)/RefPar.MatProps.tau*100;

%% Calcualte error statistics
%Bulk Modulus
sysErr.G=squeeze(mean(Errors.G,3));
ranErr.G=squeeze(std(Errors.G,0,3));
totErr.G=abs(sysErr.G)+2*ranErr.G;



% Time Constant
sysErr.tau=squeeze(mean(Errors.tau,3));
ranErr.tau=squeeze(std(Errors.tau,0,3));
totErr.tau=abs(sysErr.tau)+2*ranErr.tau;

%% Find location with the minimum shear modulus identification error
Ident.Gmean=squeeze(mean(Ident.G,3));

Ident.Gind=find(totErr.G==min(totErr.G,[],'all'));
[Ident.GtempInd,Ident.GspaInd]=ind2sub(size(totErr.G),Ident.Gind);
Ident.GminErr=Ident.Gmean(totErr.G==min(totErr.G,[],'all'));
Ident.OptSpaKern=SpaKernVec(Ident.GspaInd);
Ident.OptTempKern=TempKernVec(Ident.GtempInd);
%% Save Results
    save(strcat(savePath,'/',TestDeg,...
        '_CombinedSpatialTemporal_Noisy','_ShearOnly.mat'),...
        'intMatrix','ub','lb','ConstitutiveParam','SolveOpts','imageNoise',...
        'CondOpts','Errors','smoothOpts','TempKernVec','SpaKernVec', ...
        'phi','extrapOpts','Ident','sysErr','ranErr','totErr');

%% Plot Systematic Results image SC
figure('units','Centimeters','outerposition',[0 0 17.8 9])

subplot(1,2,1)
imagesc(SpaKernVec,TempKernVec,sysErr.G)
xlabel('Spatial Kernel (pixels)')
ylabel('Temporal Kernel (frames)')
title('G_1 Error_{sys} Identified G (%)')
crameri('-roma')
colorbar
set(gca,'YDir','normal')
caxis([-25,25])

subplot(1,2,2)
imagesc(SpaKernVec,TempKernVec,sysErr.tau)
xlabel('Spatial Kernel (pixels)')
ylabel('Temporal Kernel (frames)')
title('\tau_1  Error_{sys} identified G (%)')
crameri('-roma')
colorbar
set(gca,'YDir','normal')
caxis([-25,25])
FigSaveName=strcat(savePath,'/',TestDeg,...
    '_NoNoise_Gonly_sysErr_imagesc');

saveas(gcf,FigSaveName,'fig')
saveas(gcf,FigSaveName,'svg')

%% Plot Results image SC
figure('units','Centimeters','outerposition',[0 0 17.8 9])
subplot(1,2,1)
imagesc(SpaKernVec,TempKernVec,sysErr.G)
xlabel('Spatial Kernel (pixels)')
ylabel('Temporal Kernel (frames)')
title('G_1 Error_{sys} Identified G (%)')
crameri('-roma')
colorbar
set(gca,'YDir','normal')
caxis([-5,5])

subplot(1,2,2)
imagesc(SpaKernVec,TempKernVec,sysErr.tau)
xlabel('Spatial Kernel (pixels)')
ylabel('Temporal Kernel (frames)')
title('\tau_1 Error_{sys} Identified G (%)')
crameri('-roma')
colorbar
set(gca,'YDir','normal')
caxis([-5,5])
FigSaveName=strcat(savePath,'/',TestDeg,...
    '_NoNoise_Gonly_sysErr_imagesc_zoomed');

saveas(gcf,FigSaveName,'fig')
saveas(gcf,FigSaveName,'svg')

%% Plot Results
figure('units','Centimeters','outerposition',[0 0 17.8 9])

subplot(1,2,1)
contourf(SpaKernVec,TempKernVec,sysErr.G,50,'LineStyle','none')
xlabel('Spatial Kernel (pixels)')
ylabel('Temporal Kernel (frames)')
title('G_1 Err_{sys} Identified G (%)')
crameri('-roma')
clim([-25,25])
colorbar

subplot(1,2,2)
contourf(SpaKernVec,TempKernVec,sysErr.tau,50,'LineStyle','none')
xlabel('Spatial Kernel (pixels)')
ylabel('Temporal Kernel (frames)')
title('\tau_1 Err_{sys} Identified G (%)')
crameri('-roma')
clim([-25,25])
colorbar

FigSaveName=strcat(savePath,'/',TestDeg,'_NoNoise_Gonly_sysErr_BWYRK');
saveas(gcf,FigSaveName,'fig')
saveas(gcf,FigSaveName,'svg')

%% Plot Random Results image SC
figure('units','Centimeters','outerposition',[0 0 17.8 9])

subplot(1,2,1)
imagesc(SpaKernVec,TempKernVec,ranErr.G)
xlabel('Spatial Kernel (pixels)')
ylabel('Temporal Kernel (frames)')
title('G_1 Error_{ran} Identified G (%)')
crameri('-roma')
colorbar
set(gca,'YDir','normal')
caxis([-25,25])

subplot(1,2,2)
imagesc(SpaKernVec,TempKernVec,ranErr.tau)
xlabel('Spatial Kernel (pixels)')
ylabel('Temporal Kernel (frames)')
title('\tau_1  Error_{ran} identified G (%)')
crameri('-roma')
colorbar
set(gca,'YDir','normal')
caxis([-25,25])
FigSaveName=strcat(savePath,'/',TestDeg,...
    '_NoNoise_Gonly_ranErr_imagesc');

saveas(gcf,FigSaveName,'fig')
saveas(gcf,FigSaveName,'svg')

%% Plot Results image SC
figure('units','Centimeters','outerposition',[0 0 17.8 9])

subplot(1,2,1)
imagesc(SpaKernVec,TempKernVec,ranErr.G)
xlabel('Spatial Kernel (pixels)')
ylabel('Temporal Kernel (frames)')
title('G_1 Error_{ran} Identified G (%)')
crameri('-roma')
colorbar
set(gca,'YDir','normal')
caxis([-5,5])

subplot(1,2,2)
imagesc(SpaKernVec,TempKernVec,ranErr.tau)
xlabel('Spatial Kernel (pixels)')
ylabel('Temporal Kernel (frames)')
title('\tau_1 Error_{ran} Identified G (%)')
crameri('-roma')
colorbar
set(gca,'YDir','normal')
caxis([-5,5])
FigSaveName=strcat(savePath,'/',TestDeg,...
    '_NoNoise_Gonly_ranErr_imagesc_zoomed');

saveas(gcf,FigSaveName,'fig')
saveas(gcf,FigSaveName,'svg')

%% Plot Results
figure('units','Centimeters','outerposition',[0 0 17.8 9])
subplot(1,2,1)
contourf(SpaKernVec,TempKernVec,ranErr.G,50,'LineStyle','none')
xlabel('Spatial Kernel (pixels)')
ylabel('Temporal Kernel (frames)')
title('G_1 Err_{ran} Identified G (%)')
crameri('-roma')
clim([-25,25])
colorbar

subplot(1,2,2)
contourf(SpaKernVec,TempKernVec,ranErr.tau,50,'LineStyle','none')
xlabel('Spatial Kernel (pixels)')
ylabel('Temporal Kernel (frames)')
title('\tau_1 Err_{ran} Identified G (%)')
crameri('-roma')
clim([-25,25])
colorbar

FigSaveName=strcat(savePath,'/',TestDeg,'_NoNoise_Gonly_ranErr_BWYRK');
saveas(gcf,FigSaveName,'fig')
saveas(gcf,FigSaveName,'svg')
%% Plot total Results image SC
figure('units','Centimeters','outerposition',[0 0 17.8 9])

subplot(1,2,1)
imagesc(SpaKernVec,TempKernVec,totErr.G)
xlabel('Spatial Kernel (pixels)')
ylabel('Temporal Kernel (frames)')
title('G_1 Error_{tot} Identified G (%)')
crameri('-roma')
colorbar
set(gca,'YDir','normal')
caxis([-25,25])

subplot(1,2,2)
imagesc(SpaKernVec,TempKernVec,totErr.tau)
xlabel('Spatial Kernel (pixels)')
ylabel('Temporal Kernel (frames)')
title('\tau_1  Error_{tot} identified G (%)')
crameri('-roma')
colorbar
set(gca,'YDir','normal')
caxis([-25,25])
FigSaveName=strcat(savePath,'/',TestDeg,...
    '_NoNoise_Gonly_totErr_imagesc');

saveas(gcf,FigSaveName,'fig')
saveas(gcf,FigSaveName,'svg')

%% Plot Results image SC
figure('units','Centimeters','outerposition',[0 0 17.8 9])

subplot(1,2,1)
imagesc(SpaKernVec,TempKernVec,totErr.G)
xlabel('Spatial Kernel (pixels)')
ylabel('Temporal Kernel (frames)')
title('G_1 Error_{tot} Identified G (%)')
crameri('-roma')
colorbar
set(gca,'YDir','normal')
caxis([-5,5])

subplot(1,2,2)
imagesc(SpaKernVec,TempKernVec,totErr.tau)
xlabel('Spatial Kernel (pixels)')
ylabel('Temporal Kernel (frames)')
title('\tau_1 Error_{tot} Identified G (%)')
crameri('-roma')
colorbar
set(gca,'YDir','normal')
caxis([-5,5])
FigSaveName=strcat(savePath,'/',TestDeg,...
    '_NoNoise_Gonly_totErr_imagesc_zoomed');

saveas(gcf,FigSaveName,'fig')
saveas(gcf,FigSaveName,'svg')

%% Plot Results
figure('units','Centimeters','outerposition',[0 0 17.8 9])

subplot(1,2,1)
contourf(SpaKernVec,TempKernVec,totErr.G,50,'LineStyle','none')
xlabel('Spatial Kernel (pixels)')
ylabel('Temporal Kernel (frames)')
title('G_1 Err_{tot} Identified G (%)')
crameri('-roma')
clim([-25,25])
colorbar

subplot(1,2,2)
contourf(SpaKernVec,TempKernVec,totErr.tau,50,'LineStyle','none')
xlabel('Spatial Kernel (pixels)')
ylabel('Temporal Kernel (frames)')
title('\tau_1 Err_{tot} Identified G (%)')
crameri('-roma')
clim([-25,25])
colorbar

FigSaveName=strcat(savePath,'/',TestDeg,'_NoNoise_Gonly_totErr_BWYRK');
saveas(gcf,FigSaveName,'fig')
saveas(gcf,FigSaveName,'svg')

%% Calculate minimum value of the cost function
MinPhi=min(phi,[],3);
% %% Load Interpolated FE data
% [IntName,IntPath]=uigetfile('', ...
%     'Choose File containing Interpolated FE data');
% Int=load(strcat(IntPath,'/',IntName));

% %% Set up variables for RMSE calculating Parfor Loop
% 
% 
% strainxx=Int.strain.x;
% strainyy=Int.strain.y;
% strainxy=Int.strain.s;
% 
% stressRef=squeeze(mean(Int.stress.x));
% timevec=timefull.vec;
% RMSEi=zeros(size(Ident.K.I));
% RMSEe=zeros(size(Ident.K.I));
% 
% tempMatE=cell(size(Ident.K.E));
% tempMatI=cell(size(Ident.K.I));
% for k=1:TempKernNum
%     for m=1:SpaKernNum
%         tempPropsE=RefPar.MatProps;
%         tempPropsE.Ki=Ident.K.E(k,m);
%         tempPropsE.tau=Ident.tau.E(k,m);
%         tempMatE{k,m}=tempPropsE;
% 
%         tempPropsI=RefPar.MatProps;
%         tempPropsI.Gi=exactPropsI.Gi;
%         tempPropsI.Ki=Ident.K.I(k,m);
%         tempPropsI.tau=Ident.tau(k,m);
%         tempMatI{k,m}=tempPropsI;
%     end
% end
% % %% Calculate RMSE Error
% % NumPts=numel(Ident.K.I);
% % fprintf('Calculating RMSE Error \n')
% % parfor n=1:NumPts
% %    itPropsI=cell2mat(tempMatI(n));
% %    itPropsE=cell2mat(tempMatE(n));
% %    
% %    TempModelE=func_ViscoConstitutiveV6(strainxx,strainyy,strainxy,...
% %         timevec,itPropsE,0,0,0);
% %     RMSEe(n)=func_calcFFRMSE(stressRef,TempModelE.Avxx);
% % 
% %    TempModelI=func_ViscoConstitutiveV6(strainxx,strainyy,strainxy,...
% %         timevec,itPropsI,0,0,0);
% %     RMSEi(n)=func_calcFFRMSE(stressRef,TempModelI.Avxx);    
% % end
% % 
% % RMSE.E=RMSEe;
% % RMSENorm.E=RMSE.E/max(abs(stressRef),[],"all");
% % 
% % RMSE.I=RMSEi;
% % RMSENorm.I=RMSE.I/max(abs(stressRef),[],"all");
% % 
% % clear RMSEe RMSEi
% % %% Plot Cost function and Normalized RMSE
% % figure('units','Centimeters','outerposition',[0 0 17.8 9])
% % subplot(2,2,1)
% % imagesc(SpaKernVec,TempKernVec,RMSENorm.E)
% % xlabel('Spatial Kernel (pixels)')
% % ylabel('Temporal Kernel (frames)')
% % title('RMSE G_{ref} (\sigma_{xx})/max(\sigma_{xx})')
% % colormap(WYRK_map)
% % cx=colorbar;
% % cx.Label.String='%';
% % set(gca,'YDir','normal')
% % 
% % 
% % subplot(2,2,2)
% % imagesc(SpaKernVec,TempKernVec,MinPhi.E)
% % xlabel('Spatial Kernel (pixels)')
% % ylabel('Temporal Kernel (frames)')
% % title('\phi_K G_{ident}')
% % colormap(WYRK_map)
% % cx=colorbar;
% % cx.Label.String='\phi';
% % 
% % set(gca,'YDir','normal')
% % 
% % subplot(1,2,1)
% % imagesc(SpaKernVec,TempKernVec,RMSENorm.I)
% % xlabel('Spatial Kernel (pixels)')
% % ylabel('Temporal Kernel (frames)')
% % title('RMSE G_{ref} (\sigma_{xx})/max(\sigma_{xx})')
% % colormap(WYRK_map)
% % cx=colorbar;
% % cx.Label.String='%';
% % set(gca,'YDir','normal')
% % 
% % 
% % subplot(1,2,2)
% % imagesc(SpaKernVec,TempKernVec,MinPhi.I)
% % xlabel('Spatial Kernel (pixels)')
% % ylabel('Temporal Kernel (frames)')
% % title('\phi_K G_{ident}')
% % colormap(WYRK_map)
% % cx=colorbar;
% % cx.Label.String='\phi';
% % 
% % set(gca,'YDir','normal')
% % 
% % FigSaveName=strcat(savePath,'/',TestDeg,...
% %     '_NoNoise_Gonly_RMSE_Phi');
% % 
% % saveas(gcf,FigSaveName,'fig')
% % saveas(gcf,FigSaveName,'svg')
% % 
% % %% Zoom in
% % figure('units','Centimeters','outerposition',[0 0 17.8 9])
% % subplot(2,2,1)
% % imagesc(SpaKernVec,TempKernVec,RMSENorm.E)
% % xlabel('Spatial Kernel (pixels)')
% % ylabel('Temporal Kernel (frames)')
% % title('RMSE G_{ref} (\sigma_{xx})/max(\sigma_{xx})')
% % colormap(WYRK_map)
% % cx=colorbar;
% % cx.Label.String='%';
% % set(gca,'YDir','normal')
% % 
% % 
% % subplot(2,2,2)
% % imagesc(SpaKernVec,TempKernVec,MinPhi.E)
% % xlabel('Spatial Kernel (pixels)')
% % ylabel('Temporal Kernel (frames)')
% % title('\phi_K G_{ident}')
% % colormap(WYRK_map)
% % cx=colorbar;
% % cx.Label.String='\phi';
% % 
% % set(gca,'YDir','normal')
% % 
% % subplot(1,2,1)
% % imagesc(SpaKernVec,TempKernVec,RMSENorm.I)
% % xlabel('Spatial Kernel (pixels)')
% % ylabel('Temporal Kernel (frames)')
% % title('RMSE G_{ident} (\sigma_{xx})/max(\sigma_{xx})')
% % colormap(WYRK_map)
% % cx=colorbar;
% % cx.Label.String='%';
% % set(gca,'YDir','normal')
% % 
% % 
% % subplot(1,2,2)
% % imagesc(SpaKernVec,TempKernVec,MinPhi.I)
% % xlabel('Spatial Kernel (pixels)')
% % ylabel('Temporal Kernel (frames)')
% % title('\phi_K G_{ident}')
% % colormap(WYRK_map)
% % cx=colorbar;
% % cx.Label.String='\phi';
% % 
% % set(gca,'YDir','normal')
% % 
% % FigSaveName=strcat(savePath,'/',TestDeg,...
% %     '_NoNoise_Gonly_RMSE_Phi_zoomed');
% % 
% % saveas(gcf,FigSaveName,'fig')
% % saveas(gcf,FigSaveName,'svg')

