function [SubStrain,SubDisp,SubAccel,...
    SubRecord] =func_SubEdgeFieldsV1(SubParam,FE,GM)
%This function is intended to replace the data, extracted via the grid
%method for in-plane displacement resolution, from syntehtically deformed
%grids with pure finite element data interpolated to the pixel coordianates 

% Author: Andrew Matejunas
%Date Completed:
%Date Verified:
%Change log:
    %V1- Changed interpolation method from self coded to interp2 native
        %function in matlab

%Function input arguments
    %SubParam- guiding parameters governing the substitution with fields:
        %fields- defines what combination of strain, displacement, and
            %acceleeration will be substituted
        %num- number of pixels from the image edges over which the
            %substitutions will be performed
        %StrainSub, AccelSub, DispSup- logical arguments defining whether
            %substitution will be perfomed on strain displacement and
            %acceleration fileds. yes if true, no if false
    %FE- structure of deformation data esxtracted directly form the finite
         %element simulation with fields
        %strain- structure of strains in x (strain.x), y (strain.y) and
            %in-plane shear (strain.s). Array sizes [ycoords,xcoords,time]
        %disp- displacement fields in x (disp.x) and y (disp.y)
        %accel- acceleration fields in x (accel.x) and y (accel.y)
        %X_vec- Vector of X coordinates
        %Y_vec- vecotr of y coordinates
    %GM- Structure containing the defomation data extrtacted from
        %synthetically deformed grid images using the grid method. Contains
        %the same fields as FE along with 
            %grid- contains information about the grid image including
                %name- the name of the particular grid (not important)
                %pitch- pitch of the grid in meters (apatial length of one
                    %period)
                %pxPerPeriod- the number of pixels in one period of the
                    %grid
                %rotAngle- Angle of rotation of the grid
                %length- Specimen length in meters
                %height- specimen height in meters
                %mPerPx- size of a pixel in meters
                %numXperiods- number of grid periods in the x direction
                %asymmPitch- 0 if x pitch and y pitch are the same, 1 if
                    %different
                %pitchX- grid pitch in the x direction
                %pitchY-grid pitch in the y direction
                
                
%Function Output arguments
    %SubStrain- strain fields with substitutions on the edges
    %SubDisp- displacement fields with substitutions on the edges
    %SubAccel- acceleration fields with substitutions on the edges
    %SubRecord- Record of the substitution parameters
    

subnum=SubParam.num;
tpts=length(GM.time.vec); %#of time data points




%% Define specimen edges
    %Note this is computationally inefficient and I can probably axe it in
        %a later version
% Opposite of impact 
LeftEdge.GM=GM.X_vec(1);
LeftEdge.FE=FE.X_vec(1);

% Impact edge
RightEdge.GM=GM.X_vec(end);
RightEdge.FE=FE.X_vec(end);

% Top Edge
TopEdge.GM=GM.Y_vec(end);
TopEdge.FE=FE.Y_vec(end);

% Bottom Edge
BEdge.GM=GM.Y_vec(1);
BEdge.FE=FE.Y_vec(1);


%% Define coordinates over which substitutions will be performed

SubCoords.B=GM.Y_vec(1:subnum); %bottom
SubCoords.T=GM.Y_vec((end-subnum+1):end); %top
SubCoords.L=GM.X_vec(1:subnum); %left (Free)
SubCoords.R=GM.X_vec((end-subnum+1):end); %righrt (impact

%% Get vectors of FE data points needed for interpolation
FEpoints.B=FE.Y_vec(FE.Y_vec<=SubCoords.B(end)); %bottom
FEpoints.T=FE.Y_vec(FE.Y_vec>=SubCoords.T(1)); %top
FEpoints.L=FE.X_vec(FE.X_vec<=SubCoords.L(end)); %left
FEpoints.R=FE.X_vec(FE.X_vec>=SubCoords.R(1)); %right

%% Identify appropriate indexes on GM data
SubIndexes.T=(length(GM.Y_vec)-subnum+1):length(GM.Y_vec);
SubIndexes.B=1:subnum;
SubIndexes.L=1:subnum;
SubIndexes.R=(length(GM.X_vec)-subnum+1):length(GM.X_vec);

%% Identify indexes that matter on FE data
FEindexes.T=(length(FE.Y_vec)-length(FEpoints.T)+1):length(FE.Y_vec);
FEindexes.B=1:length(FEpoints.B);
FEindexes.L=1:length(FEpoints.L);
FEindexes.R=(length(FE.X_vec)-length(FEpoints.R)+1):length(FE.X_vec);


%% initialize variables

%displacements
SubDispX=GM.disp.x;
SubDispY=GM.disp.y;
RDispX=FE.disp.x;
RDispY=FE.disp.y;

SubStrainXX=GM.strain.x;
SubStrainYY=GM.strain.y;
SubStrainXY=GM.strain.s;
RStrainXX=FE.strain.x;
RStrainYY=FE.strain.y;
RStrainXY=FE.strain.s;

SubAccelX=GM.accel.x;
SubAccelY=GM.accel.y;
RAccelX=FE.accel.x;
RAccelY=FE.accel.y;



%% Perform interpolaton on top data points
%substitute directly on top Edge

%Define interpolation grid points
X=FEpoints.T; %Ycoordinates of FE referance data
Y=FE.X_vec;   %X coordinates of FE reference data
Yq=GM.X_vec;

for k=1:tpts
    %% Displacement 
    if SubParam.DispSub==true
        %X direction
     VX=squeeze(RDispX(FEindexes.T,:,k));
    SubDispX(SubIndexes.T,:,k)=interp2(Y,X,VX,SubCoords.T',Yq)';
    %Ydirection
    VY=squeeze(RDispY(FEindexes.T,:,k));
    SubDispY(SubIndexes.T,:,k)=interp2(Y,X,VY,SubCoords.T',Yq)';
    end
    
    %% Strain
    if SubParam.StrainSub==true
        %XX
    VXX=squeeze(RStrainXX(FEindexes.T,:,k));
    SubStrainXX(SubIndexes.T,:,k)=interp2(Y,X,VXX,SubCoords.T',Yq)';
        %YY
    VYY=squeeze(RStrainYY(FEindexes.T,:,k));
    SubStrainYY(SubIndexes.T,:,k)=interp2(Y,X,VYY,SubCoords.T',Yq)';
        %Shear
    VXY=squeeze(RStrainXY(FEindexes.T,:,k));
    SubStrainXY(SubIndexes.T,:,k)=interp2(Y,X,VXY,SubCoords.T',Yq)';
    end
    
    %% Acceleration
     if SubParam.AccelSub==true
        %X direction
    VX=squeeze(RAccelX(FEindexes.T,:,k));
    SubAccelX(SubIndexes.T,:,k)=interp2(Y,X,VX,SubCoords.T',Yq)';
        %Ydirection
    VY=squeeze(RAccelY(FEindexes.T,:,k));
    SubAccelY(SubIndexes.T,:,k)=interp2(Y,X,VY,SubCoords.T',Yq)';
     end
end

%% Perform interpolaton on Bottom data points
%substitute directly on top Edge

%Define interpolation grid points
X=FEpoints.B; %Ycoordinates of FE referance data
Y=FE.X_vec;   %X coordinates of FE reference data
Yq=GM.X_vec;

for k=1:tpts
    %% Displacement 
    if SubParam.DispSub==true
        %X direction
     VX=squeeze(RDispX(FEindexes.B,:,k));
    SubDispX(SubIndexes.B,:,k)=interp2(Y,X,VX,SubCoords.B',Yq)';
    %Ydirection
    VY=squeeze(RDispY(FEindexes.B,:,k));
    SubDispY(SubIndexes.B,:,k)=interp2(Y,X,VY,SubCoords.B',Yq)';
    end
    
    %% Strain
    if SubParam.StrainSub==true
        %XX
    VXX=squeeze(RStrainXX(FEindexes.B,:,k));
    SubStrainXX(SubIndexes.B,:,k)=interp2(Y,X,VXX,SubCoords.B',Yq)';
        %YY
    VYY=squeeze(RStrainYY(FEindexes.B,:,k));
    SubStrainYY(SubIndexes.B,:,k)=interp2(Y,X,VYY,SubCoords.B',Yq)';
        %Shear
    VXY=squeeze(RStrainXY(FEindexes.B,:,k));
    SubStrainXY(SubIndexes.B,:,k)=interp2(Y,X,VXY,SubCoords.B',Yq)';
    end
    
    %% Acceleration
     if SubParam.AccelSub==true
        %X direction
    VX=squeeze(RAccelX(FEindexes.B,:,k));
    SubAccelX(SubIndexes.B,:,k)=interp2(Y,X,VX,SubCoords.B',Yq)';
        %Ydirection
    VY=squeeze(RAccelY(FEindexes.B,:,k));
    SubAccelY(SubIndexes.B,:,k)=interp2(Y,X,VY,SubCoords.B',Yq)';
     end
end

%% Perform interpolaton on Left data points
%substitute directly on top Edge

%Define interpolation grid points
X=FE.Y_vec; %Ycoordinates of FE referance data
Y=FEpoints.L;   %X coordinates of FE reference data
Xq=GM.Y_vec;

for k=1:tpts
    %% Displacement 
    if SubParam.DispSub==true
        %X direction
     VX=squeeze(RDispX(:,FEindexes.L,k));
    SubDispX(:,SubIndexes.L,k)=interp2(Y,X,VX,Xq,SubCoords.L')';
    %Ydirection
    VY=squeeze(RDispY(:,FEindexes.L,k));
    SubDispY(:,SubIndexes.L,k)=interp2(Y,X,VY,Xq,SubCoords.L')';
    end
    
    %% Strain
    if SubParam.StrainSub==true
        %XX
    VXX=squeeze(RStrainXX(:,FEindexes.L,k));
    SubStrainXX(:,SubIndexes.L,k)=interp2(Y,X,VXX,Xq,SubCoords.L')';
        %YY
    VYY=squeeze(RStrainYY(:,FEindexes.L,k));
    SubStrainYY(:,SubIndexes.L,k)=interp2(Y,X,VYY,Xq,SubCoords.L')';
        %Shear
    VXY=squeeze(RStrainXY(:,FEindexes.L,k));
    SubStrainXY(:,SubIndexes.L,k)=interp2(Y,X,VXY,Xq,SubCoords.L')';
    end
    
    %% Acceleration
     if SubParam.AccelSub==true
        %X direction
    VX=squeeze(RAccelX(:,FEindexes.L,k));
    SubAccelX(:,SubIndexes.L,k)=interp2(Y,X,VX,Xq,SubCoords.L')';
        %Ydirection
    VY=squeeze(RAccelY(:,FEindexes.L,k));
    SubAccelY(:,SubIndexes.L,k)=interp2(Y,X,VY,Xq,SubCoords.L')';
     end
end

%% Perform interpolaton on Right data points
%substitute directly on top Edge

%Define interpolation grid points
X=FE.Y_vec; %Ycoordinates of FE referance data
Y=FEpoints.R;   %X coordinates of FE reference data
Xq=GM.Y_vec;

for k=1:tpts
    %% Displacement 
    if SubParam.DispSub==true
        %X direction
     VX=squeeze(RDispX(:,FEindexes.R,k));
    SubDispX(:,SubIndexes.R,k)=interp2(Y,X,VX,Xq,SubCoords.R')';
    %Ydirection
    VY=squeeze(RDispY(:,FEindexes.R,k));
    SubDispY(:,SubIndexes.R,k)=interp2(Y,X,VY,Xq,SubCoords.R')';
    end
    
    %% Strain
    if SubParam.StrainSub==true
        %XX
    VXX=squeeze(RStrainXX(:,FEindexes.R,k));
    SubStrainXX(:,SubIndexes.R,k)=interp2(Y,X,VXX,Xq,SubCoords.R')';
        %YY
    VYY=squeeze(RStrainYY(:,FEindexes.R,k));
    SubStrainYY(:,SubIndexes.R,k)=interp2(Y,X,VYY,Xq,SubCoords.R')';
        %Shear
    VXY=squeeze(RStrainXY(:,FEindexes.R,k));
    SubStrainXY(:,SubIndexes.R,k)=interp2(Y,X,VXY,Xq,SubCoords.R')';
    end
    
    %% Acceleration
     if SubParam.AccelSub==true
        %X direction
    VX=squeeze(RAccelX(:,FEindexes.R,k));
    SubAccelX(:,SubIndexes.R,k)=interp2(Y,X,VX,Xq,SubCoords.R')';
        %Ydirection
    VY=squeeze(RAccelY(:,FEindexes.R,k));
    SubAccelY(:,SubIndexes.R,k)=interp2(Y,X,VY,Xq,SubCoords.R')';
     end
end

%% Assign function outputs
SubRecord.Param=SubParam;
SubRecord.Coords=SubCoords;
SubRecord.GMindexes=SubIndexes;
SubRecord.FEindexes=FEindexes;

SubStrain.x=SubStrainXX;
SubStrain.y=SubStrainYY;
SubStrain.s=SubStrainXY;

SubDisp.x=SubDispX;
SubDisp.y=SubDispY;

SubAccel.x=SubAccelX;
SubAccel.y=SubAccelY;

end

