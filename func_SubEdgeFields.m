function [SubStrain,SubDisp,SubAccel,...
    SubRecord] =func_SubEdgeFields(SubParam,FE,GM)
%This function is intended to replace the data, extracted via the grid
%method for in-plane displacement resolution, from syntehtically deformed
%grids with pure finite element data interpolated to the pixel coordianates 

% Author: Andrew Matejunas
%Date Completed:
%Date Verified:
%Change log:

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

%% Perform interpolaton on top data points
%initialize variables

SubDispX=GM.disp.x;
SubDispY=GM.disp.y;

SubStrainXX=GM.strain.x;
SubStrainYY=GM.strain.y;
SubStrainXY=GM.strain.s;

SubAccelX=GM.accel.x;
SubAccelY=GM.accel.y;

%substitute directly on top Edge

for k=1:length(SubIndexes.T)-1
    %% Displacement 
    if SubParam.DispSub==true
    
    SubInd=SubIndexes.T(k);
    
    FEgreater=FEpoints.T(FEpoints.T>=SubCoords.T(k));
    FEgInd=FEindexes.T(FEpoints.T>=SubCoords.T(k));
    
    FEless=FEpoints.T(FEpoints<=SubCoords.T(k));
    FElind=FEindexes.T(FEpoints.T<=SubCoords.T(k));
    
    end    
end

SubRecord.Param=SubParam;
Subrecord.Coords=SubCoords;
Subrecord.GMindexes=SubIndexes;
Subrecord.FEindexes=FEindexes;



end

