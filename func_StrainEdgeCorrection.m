function [disp,strain,accel]=func_StrainEdgeCorrection(grid,pos,...
    disp,strain,accel,time,...
    corOpts)

% This script is written to correct for errors caused by grid method
    % processing that accor along the edges of the grid. The method
    % linearly interpolates the kinematic field between the free surface
    % and the closest good datapoint. 
    
 % Author: Andrew Matejunasco
 %date completed:
 %date verified:
 
 %Change log:
 
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
 
 % Function output arguments
    %disp- corrected displacement field
    %strain- corrected strain fields
    %accel- correscted acceleration fields
 pitch=grid.pxPerPeriod;    
 Xcoords=pos.x;
 Ycoords=pos.y;
 CorNum=corOpts.CorNum;
 
% correct strain only


    %% set shear and YY strains to 0 on top and bottom free edges
    %shear
    strain.s(1,:,:)=0;
    strain.s(end,:,:)=0;
    
    %YY
    strain.y(1,:,:)=0;
    strain.y(end,:,:)=0;
    
 %% Perform linear interpolation
  %% Record coordinates of interest
  %coordinate of bottom free edge
  bottomcoord=Ycoords(1);
  %first valid coordinate from the bottom free edge
  Bvalidcoord=Ycoords(1+CorNum*pitch);
  
  %coordinate of top free edge and first valid point
  topcoord=Ycoords(end);
  Tvalidcoord=Ycoords(end-CorNum*pitch);
  
  %Coordinates to be interpolated
  Binterp=Ycoords(2:CorNum*pitch);
  Tinterp=Ycoords((end-1):(end-CorNum*pitch+1));
 %% record strains for bottom free edge 
   BottomstrainXY=strain.s(1,:,:);
   BottomstrainYY=strain.y(1,:,:);
   
   BvalidstrainXY=strain.s(1+CorNum*pitch,:,:);
   BvalidstrainYY=strain.y(1+CorNum*pitch,:,:);
   
  
   %% record strains for top free edge
   TopstrainXY=strain.s(end,:,:);
   TopstrainYY=strain.y(end,:,:);
    
   TvalidstrainXY=strain.s(end-CorNum*pitch,:,:);
   TvalidstrainYY=strain.y(end-CorNum*pitch,:,:);
  
  %% Perform the interpolation
  for n=1:length(time.vec)
      for m=1:length(Xcoords)
  %% Bottom Free edge
 strain.s(2:(CorNum*pitch),m,n)=interp1([bottomcoord,Bvalidcoord],...
         [BottomstrainXY(1,m,n),BvalidstrainXY(1,m,n)],...
         Binterp);
 strain.y(2:CorNum*pitch,m,n)=interp1([bottomcoord,Bvalidcoord],...
         [BottomstrainYY(1,m,n),BvalidstrainYY(1,m,n)],...
         Binterp);    
  %% Top Free edge
 strain.s((end-1):(end-CorNum*pitch+1),m,n)=interp1([topcoord,Tvalidcoord],...
         [TopstrainXY(1,m,n),TvalidstrainXY(1,m,n)],...
         Tinterp);
strain.y((end-1):(end-CorNum*pitch+1),m,n)=interp1([topcoord,Tvalidcoord],...
         [TopstrainYY(1,m,n),TvalidstrainYY(1,m,n)],...
         Tinterp);
      end
  end
  
   
   

     
    
 end
     
