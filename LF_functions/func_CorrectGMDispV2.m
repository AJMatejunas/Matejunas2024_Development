function[disp,DispCorr,grid,ProgramVersions]=func_CorrectGMDispV2(disp,DispCorr,...
    grid,pos)
%% Image Based Inertial Impact (IBII) Test Grid Method Displacement
    %% Correction

%This function is intended to correct displacement fields output by the
%Grid Method algorith. The grid method method produces erroneous results
%within one grid pitch of the free edges. Additionally, rotations induced
%by shear deformation can propagate error spikes well beyond a single grid
%pitch into the specimen, which also need to be corrected.
    %Note: No corrections are performed on the impact edge as that data can
    %be fully excluded from the VFM algorithm

%Author: Andrew Matejunas

%Date Created:
%Date Verified:

%Change log/ version history

%define program version used for correction
ProgramVersions.Correction_Algorithm='CorrectGMDisp';

%Function input arguments
    %grid- structure defining grid parameters. Input from the grid method
        %processing algorithm
    %disp- Raw displacement fields output by the grid method processing
        %algorithm NOTE: this version of the algorithm performs correction
        %directly on the displacement fields prior to calculation of
        %acceleration and strain fields. Thus, the strain and acceleration
        %fields are not input into this algorithm.
            %Matrix size- [Ycoordinates,Xcoordinates,Timeintervals]
    %DispCorr- Structure containing the parameters for performing
        %displacement edge corrections with fields
            %Opt- Determines whether correction will be performed
            %Method- Determines how correction will be performed with
                %options
                    %Direct- Sets displacement at edges equal to last valid
                        %value (likely not a good approximation
                    %LinGrad- linearly extrapolates dispalcement from the
                        %last valid value (may need to add higher order
                        %extrapolations)
             %int- number of pixels that will be corrected
             %PitchFitKern=number of grid pitches to obtain fit slope
     %pos-structure containing spatial data with fields
        %specimenLoc- location of the bottom left and top right corners of
            %a rectangular specimen
        %x- x coordinates in vector form
        %y- y coordinates in vector form
        %xGrid- X coordinate of every point (matrix form)
        %yGrid- Y coordinate of every point (matrix form)
        %xStep- interval between x coordinate in m
        %yStep- interval between y coordinates in m
                        
%Function output arguments
    %disp- Corrected Displacement fields
    %DispCorr- Stucture containing parameters for performing displacement
        %corrections
            %DispCorr.CorNum- added field for number of grid pitches over
                %which the correction has been performed
                %ggrid- structure describing grid method parameters
    %ProgramVersions- Data structure for tracking which programs and
        %algorthms were used in data processing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Temporarily save raw displacementfieds in case reference is needed
Rawdisp=disp;

%% Set up needed variables for correction

%pixels in one grid pitch
pitch=grid.pxPerPeriod;

%X and y Coordinates
Xcoords=pos.x;
Ycoords=pos.y;

% Number of pts to perform displacement correction on
CorrPts=DispCorr.int;

% Save number of grid pitches that are corrected
DispCorr.CorNum=CorrPts/pitch;

%% Identify indexes of free surfaces

%size of displacement matrix
DispSize=size(disp.x);
%top index
TopIn=DispSize(1);
Tcoord=Ycoords(TopIn);

%bottom Index
Bind=1;
Bcoord=Ycoords(Bind);

%Index of left free edge (Note assumes impact is on the right edge)
Lin=1;
Lcoord=Xcoords(Lin);

%% identify indexes and coordinates of last valid data point

%Top
TVin=TopIn-CorrPts-1;
TVcoord=Ycoords(TVin);

%Bottom
BVin=Bind+CorrPts+1;
BVcoord=Ycoords(BVin);

%left
LVin=Lin+CorrPts+1;
LVcoord=Xcoords(Lin);

%% Extract last valid displacements

%Matrix of last valid displacements along top edge
TVDisp.x=disp.x(TVin,:,:); 
TVDisp.y=disp.y(TVin,:,:);

%bottom edge
BVDisp.x=disp.x(BVin,:,:);
BVDisp.y=disp.y(BVin,:,:);

%left edge
LVDisp.x=disp.x(:,LVin,:);
LVDisp.y=disp.y(:,LVin,:);

%% Perform the corrections

%% If using the direct method
switch DispCorr.Method
    case 'Direct'
        ProgramVersions.Correction_Method='Direct';
        %% correct displacements on top free surface
        for k=TVin:TopIn
            disp.x(k,:,:)=TVDisp.x;
            disp.y(k,:,:)=TVDisp.y;
        end
        %% correct on bottom free surface
        for k=Bind:BVin
            disp.x(k,:,:)=BVDisp.x;
            disp.y(k,:,:)=BVDisp.y;
        end
        %% Correct left free surface
        for k=Lin:LVin
            disp.x(:,k,:)=LVDisp.x;
            disp.y(:,k,:)=LVDisp.y;
        end
        
    case 'LinGrad'
        fprintf('Linear displacement corrections not yet implemented \n')
    %% Identify indexes and of data for linear fit
    
        fitKern=DispCorr.PitchFitKern*pitch;
        DispCorr.fitKernal=fitKern;
        
        %top Edge
        TfitIn=(TVin-fitKern+1):TVin;
        TfitDisp.x=disp.x(TfitIn,:,:);
        TfitDisp.y=disp.y(TfitIn,:,:);
        
        %Bottom Edge
        BfitIn=BVin:(BVin+fitKern-1);
        BfitDisp.x=disp.x(BfitIn,:,:);
        BfitDisp.y=disp.y(BfitIn,:,:);
        
        %Left Edge
        LfitIn=LVin:(LVin+fitKern-1);
        LfitDisp.x=disp.x(:,LfitIn,:);
        LfitDisp.y=disp.y(:,LfitIn,:);
        
       
    %% Identify Coordinates of previous (reference) period
        %Top Edge Reference 
        TfitCoord=Ycoords(TfitIn);

        %Bottom Referecne Edge
        BfitCoord=Ycoords(BfitIn);

        %Left Reference edge
        LfitCoord=Xcoords(LfitIn);

      

        
    %% Calculate Average slope of previous period
        %Top
        ItCount=0;
        for k=length(TfitIn):-1:2
            ItCount=ItCount+1;
            ind1=TfitIn(k);     %Inner
            ind2=TfitIn(k-1);     %outer
            Disp1=disp.y(ind1,:,:);  
            Disp2=disp.y(ind2,:,:);
            y1=TfitCoord(k);
            y2=TfitCoord(k-1);
            FitSlope.T(ItCount,:,:)=(Disp2-Disp1)/(y2-y1);
        end
        %calculate average slope over fit interval
        FitSlope.Tav=mean(FitSlope.T,1);
        %Bottom
        ItCount=0;
         for k=2:length(BfitIn)
             ItCount=ItCount+1;
            ind1=BfitIn(k);     %Inner
            ind2=BfitIn(k-1);     %outer
            Disp1=disp.y(ind1,:,:);  
            Disp2=disp.y(ind2,:,:);
            y1=BfitCoord(k);
            y2=BfitCoord(k-1);
            FitSlope.B(ItCount,:,:)=(Disp2-Disp1)/(y2-y1);
         end
         FitSlope.Bav=mean(FitSlope.B,1);
         %Left
         ItCount=0;
        for k=2:length(LfitIn)
            ItCount=ItCount+1;
            ind1=LfitIn(k);     %Inner
            ind2=LfitIn(k-1);     %outer
            Disp1=disp.x(:,ind1,:);  
            Disp2=disp.x(:,ind2,:);
            y1=LfitCoord(k);
            y2=LfitCoord(k-1);
            FitSlope.L(:,ItCount,:)=(Disp2-Disp1)/(y2-y1);
        end
        FitSlope.Lav=mean(FitSlope.L,2);
    %% Perfom Linear Correction
       for k=(TVin+1):TopIn
           disp.y(k,:,:)=disp.y(k-1,:,:)+FitSlope.Tav.*(Ycoords(k)...
               -Ycoords(k-1));
       end
       %Bottom
       for k=(BVin-1):-1:Bind
         disp.y(k,:,:)=disp.y(k+1,:,:)-FitSlope.Bav.*(Ycoords(k+1)...
               -Ycoords(k));
       end

       % left
       for k=(LVin-1):-1:Lin
         disp.x(:,k,:)=disp.x(:,k+1,:)-FitSlope.Lav.*(Xcoords(k+1)...
               -Xcoords(k));
       end
end
end