function [disp,strain,accel,ProgramVersions]=func_StrainEdgeCorrectionV2(grid,pos,...
    disp,strain,accel,time,...
    corOpts,ProgramVersions)
% This script is written to correct for errors caused by grid method
    % processing that occur along the edges of the grid. The method
    % linearly interpolates the kinematic field between the free surface
    % and the closest good datapoint. 
    
 % Author: Andrew Matejunas
 %date completed: 4/26/2021
 %date verified:  4/26/2021
 
 %Change log:
    % V1 4/26/2021- Hardcoded interpolation algorithm myself
    % V2 2022-02-09- Adding Displacement and strain corrections to top and
                    % bottom free edge
    
ProgramVersions.Correction_Algorithm='StrainEdgeCorrectionV2';    
 
 % Function input arguments
    %grid- structure defining grid parameters. Output by the grid
        %processing algorithm
    %strain- strain fields output by the grid processing algorithm
    %disp-   displacement fields output by the grid processing algorithm
    %accel-  acceleration fields output by the grid processing algorithm
    %time- temporal data
    %corOpts- Structure containing true/false correction options contains
        %fields
            %corDisp- if true displacement field is corrected. if false no
                %correction is perfor
            %corStrain- if true strain is corrected. If false strain is
                %recalculated from spatial differentiation of the displacement
                %field
            %CorAccel-if true acceleration field is corrected indepently of
                %displacement field. If false acceleration field is
                %recalculated from temporal differentiation
            %CorNum- number of grid pitches that will be corrected at each
            %edge
            %XXMethod- Method for correcting strain on the top
                          %and bottom free surface
                            %Direct- Assumes strain at surface=Last valid
                                %strain
                            %linear- extends the average linear gradient
                                %over the last valid period to the top and
                                %bottom free surfaces. Will not allow a
                                %sign chang
 
 % Function output arguments
      %strain- corrected strain fields
    
 pitch=grid.pxPerPeriod;    
 Xcoords=pos.x;
 Ycoords=pos.y;
 
 CorNum=corOpts.CorNum; %Number of grid pitches corrected
 CorPts=CorNum*pitch;
 
 
% correct strain only
    
    %% Identify top index
        StrainSize=size(strain.s);
        TopIn=StrainSize(1);
        
       
    %% set shear and YY strains to 0 on top and bottom free edges
    %shear
    strain.s(1,:,:)=0;
    strain.s(TopIn,:,:)=0;
    
    %YY
    strain.y(1,:,:)=0;
    strain.y(TopIn,:,:)=0;
    
 %% Perform linear interpolation on strainyy, and strain XY
  %% Record coordinates of interest
  %coordinate of bottom free edge
  bottomcoord=Ycoords(1);
  %first valid coordinate from the bottom free edge
  BVin=1+CorPts; %index of bottom valid point
  Bvalidcoord=Ycoords(BVin);
  
  %coordinate of top free edge and first valid point
  topcoord=Ycoords(TopIn);
  
  TVin=TopIn-CorPts-1;
  Tvalidcoord=Ycoords(TVin);
  
  %Coordinates to be interpolated
  Bind=1:(BVin-1);
  Binterp=Ycoords(Bind);
  Tind=(TVin+1):(TopIn-1);
  Tinterp=Ycoords(Tind);
 %% record strains for bottom free edge 
   BottomstrainXY=strain.s(1,:,:);
   BottomstrainYY=strain.y(1,:,:);
   
   
   BvalidstrainXY=strain.s(BVin,:,:);
   BvalidstrainYY=strain.y(BVin,:,:);
   
  
   %% record strains for top free edge
   TopstrainXY=strain.s(TopIn,:,:);
   TopstrainYY=strain.y(TopIn,:,:);
   
   
   TvalidstrainXY=strain.s(TVin,:,:);
   TvalidstrainYY=strain.y(TVin,:,:);
  
   %% Find slopes
   BslopeXY=BvalidstrainXY/(Ycoords(BVin)-Ycoords(1));
   BslopeYY=BvalidstrainYY/(Ycoords(BVin)-Ycoords(1));
    
   TslopeXY=TvalidstrainXY/(Ycoords(TVin)-Ycoords(TopIn));
   TslopeYY=TvalidstrainYY/(Ycoords(TVin)-Ycoords(TopIn));
   
  %% Perform the interpolation
  for n=1:length(time.vec)
      for m=1:length(Xcoords)
          for k=1:length(Bind)
  %% Bottom Free edge
 strain.s(Bind(k),m,n)=BslopeXY(1,m,n)*(Ycoords(Bind(k))-Ycoords(1));
 strain.y(Bind(k),m,n)=BslopeYY(1,m,n)*(Ycoords(Bind(k))-Ycoords(1));
          end
          
  %% Top Free edge
         for k=1:length(Tind)
 strain.s(Tind(k),m,n)=TslopeXY(1,m,n)*(Ycoords(Tind(k))-Ycoords(TopIn));
 strain.y(Tind(k),m,n)=TslopeYY(1,m,n)*(Ycoords(Tind(k))-Ycoords(TopIn));
         end 
        end
  end
  
%% Perform strainXX corection
TvalidstrainXX=strain.x(TVin,:,:); 
BvalidstrainXX=strain.s(BVin,:,:);

% if Strainxx_new=strainXX_valid
    
if corOpts.XXMethod=='Direct'
   for k=1:length(Tind)
       strain.x(Tind(k),:,:)=TvalidstrainXX;
   end
   for k=1:length(Bind)
       strain.x(Bind(k),:,:)=BvalidstrainXX;
   end
end

     
    
 end
     
