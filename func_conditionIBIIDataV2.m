function [SG,accel,strain,X_vec,time]=func_conditionIBIIDataV2(SG,accel,...
    strain,X_vec,time,CondOpts)
% This script is written to perform conditioning operations on the data
% structures used for viscoelastic constitutive parameter identification
% using the image based inertial impact test. Including Downsampling,
% censoring data, and spatial and temporal smoothing.

%Author: Andrew Matejunas

%Date Completed: 2022-06-02

%Version History/Change log:
    %2022-06-02- Original version (does not include smoothing options)
    %V2-2023/02/22- Added ability to remove data from some of the last
                    %frames in the dataset

%Function input arguments
    %SG- Structure of average stresses computed with the stress gage
        %equations with fields
            %x- Normal Stresses in the X-direction (parallel to impact)
            %s- In plane shear stresses
    %accel- Structure of accelerations with fields
            %x- x direction (parallel to impact)
            %y- y direction (perpendicular to impact)
    %strain- Structure of in plane strain data with fields
            %x- normal strains in x direction
            %y- normal strains in y direction
            %s- in plane shear strains
    %X_vec- Vector of specimen X coordinates
    %CondOpts- Conditioning optinons that determine how the other inputs
        %will be manipulated. With fieds
            %ImpCens- number of data points to be removed from impact edge
            %FreeCens- number of points to be removed from free edge
            %Xds- downsampling factor for the X data
            %Yds- downsampling factor for the Y data
            %(Smoothing options to be added in future)
            %CutEndFrames- determines what frame to remove from the end of the
                      %data set
    %Structure of time data

%Function output arguments
    %SG- Conditioned stress gage stresses 
    %accel- Conditioned accelerations
    %strain- Conditioned strain fields
    %X_vec- conditioned X coordinate vector
        
        
        
    
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% record total frames
TotFrames=length(time.vec);

%% Extract the conditioning options
ImpCens=CondOpts.ImpCens;
FreeCens=CondOpts.FreeCens;
Xds=CondOpts.Xds;
Yds=CondOpts.Yds;


%% Censor Impact End
%X_vector
X_vec((end-ImpCens):end)=[];

%censor strain
strain.x(:,(end-ImpCens):end,:)=[];
strain.y(:,(end-ImpCens):end,:)=[];
strain.s(:,(end-ImpCens):end,:)=[];

%censor accel
accel.x(:,(end-ImpCens):end,:)=[];
accel.y(:,(end-ImpCens):end,:)=[];

%censor SG
SG.x((end-ImpCens):end,:)=[];
SG.s((end-ImpCens):end,:)=[];

%% Free end
%X_vector
X_vec(1:FreeCens)=[];

%censor strain
strain.x(:,1:FreeCens,:)=[];
strain.y(:,1:FreeCens,:)=[];
strain.s(:,1:FreeCens,:)=[];

%censor accel
accel.x(:,1:FreeCens,:)=[];
accel.y(:,1:FreeCens,:)=[];

%censor SG
SG.x(1:FreeCens,:)=[];
SG.s(1:FreeCens,:)=[];

%% Downsample in the X direction
[SG,accel,strain,X_vec]=func_KinFieldsDownsampling(SG,accel,...
    strain,X_vec,time,Xds,Yds);

%% Remove frames from the end
% Save number of frames to remove
EndCens=CondOpts.CutEndFrames;
EndInt=(TotFrames-EndCens+1):TotFrames;

%Censor time vector
time.vec(EndInt)=[];

%censor strains
strain.x(:,:,EndInt)=[];
strain.y(:,:,EndInt)=[];
strain.s(:,:,EndInt)=[];

%censor accelerations
accel.x(:,:,EndInt)=[];
accel.y(:,:,EndInt)=[];

%censor stress guage data
SG.x(:,EndInt)=[];
SG.s(:,EndInt)=[];

end