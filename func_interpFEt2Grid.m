function [posInt,dispInt,accelInt,strainInt,stressInt] =...
    func_interpFEt2Grid(FEpos,GridPos,disp,accel,strain,stress,interpMethod)

%This script is written to linearlly interpolate finite element kinematic
    %fields to the coordinate system of the grid method images.

%Author: Andrew Matejunas

%Date created: 2023/03/06

%Version history/Change log
    %2023/03/08- added an input for the interpolation method

%function input arguments
    %FEpos- pos data structure from the finite element output
    %GridPos- pos data structure from grid method pocessing
    %disp- displacement structure output by the finite element simulation
    %accel- acceleration data structure output by the finite element
        %simulation
    %strain- strain data structure output by the finite element simulation
    %stress- stress data structure output by the finite element simulation

%function output arguments
    %posInt- pos data structure that FE data is interpolated over
    %dispInt- interpolated displacement data structure
    %accelInt- interpolated acceleration data structure
    %strainInt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save GridPos as the interpolated
posInt=GridPos;

%% time length
ltime=size(disp.x,3);
%% Initialize interpolated variables
[dispInt.x,dispInt.y,...
    accelInt.x,accelInt.y,...
    strainInt.x,strainInt.y,strainInt.s,...
    stressInt.x,stressInt.y,stressInt.s]=deal(...
    zeros([size(GridPos.xGrid),ltime]));
%% perform interpolation
for t=1:ltime
    %% x displacement
    xdisp=squeeze(disp.x(:,:,t));
    dispInt.x(:,:,t)=interp2(FEpos.xGrid,FEpos.yGrid,xdisp,GridPos.xGrid, ...
        GridPos.yGrid,interpMethod);

     %% y displacement
    ydisp=squeeze(disp.y(:,:,t));
    dispInt.y(:,:,t)=interp2(FEpos.xGrid,FEpos.yGrid,...
        ydisp, ...
        GridPos.xGrid,GridPos.yGrid,interpMethod);
    %% x acceleration
    xaccel=squeeze(accel.x(:,:,t));
    accelInt.x(:,:,t)=interp2(FEpos.xGrid,FEpos.yGrid,...
        xaccel, ...
        GridPos.xGrid,GridPos.yGrid,interpMethod);
    
  
    %% y acceleration
    yaccel=squeeze(accel.y(:,:,t));
    accelInt.y(:,:,t)=interp2(FEpos.xGrid,FEpos.yGrid,...
        yaccel, ...
        GridPos.xGrid,GridPos.yGrid,interpMethod);

    %% x strain
    xstrain=squeeze(strain.x(:,:,t));
    strainInt.x(:,:,t)=interp2(FEpos.xGrid,FEpos.yGrid,...
        xstrain, ...
        GridPos.xGrid,GridPos.yGrid,interpMethod);
    
  
    %% y strain
    ystrain=squeeze(strain.y(:,:,t));
    strainInt.y(:,:,t)=interp2(FEpos.xGrid,FEpos.yGrid,...
        ystrain, ...
        GridPos.xGrid,GridPos.yGrid,interpMethod);

    %% shear strain
    sstrain=squeeze(strain.s(:,:,t));
    strainInt.s(:,:,t)=interp2(FEpos.xGrid,FEpos.yGrid,...
        sstrain, ...
        GridPos.xGrid,GridPos.yGrid,interpMethod);


    %% x stress
    xstress=squeeze(stress.x(:,:,t));
    stressInt.x(:,:,t)=interp2(FEpos.xGrid,FEpos.yGrid,...
        xstress, ...
        GridPos.xGrid,GridPos.yGrid,interpMethod);
    
  
    %% y stress
    ystress=squeeze(stress.y(:,:,t));
    stressInt.y(:,:,t)=interp2(FEpos.xGrid,FEpos.yGrid,...
        ystress, ...
        GridPos.xGrid,GridPos.yGrid,interpMethod);

    %% shear strain
    sstress=squeeze(stress.s(:,:,t));
    stressInt.s(:,:,t)=interp2(FEpos.xGrid,FEpos.yGrid,...
        sstress, ...
        GridPos.xGrid,GridPos.yGrid,interpMethod);
end
%% fix NaNs
extrapOpts.fixNaNKernal=10;

dispInt.x=func_fixNaNsInDisp(extrapOpts,dispInt.x);
dispInt.y=func_fixNaNsInDisp(extrapOpts,dispInt.y);
accelInt.x=func_fixNaNsInDisp(extrapOpts,accelInt.x);
accelInt.y=func_fixNaNsInDisp(extrapOpts,accelInt.y);
strainInt.x=func_fixNaNsInDisp(extrapOpts,strainInt.x);
strainInt.y=func_fixNaNsInDisp(extrapOpts,strainInt.y);
strainInt.s=func_fixNaNsInDisp(extrapOpts,strainInt.s);
stressInt.x=func_fixNaNsInDisp(extrapOpts,stressInt.x);
stressInt.y=func_fixNaNsInDisp(extrapOpts,stressInt.y);
stressInt.s=func_fixNaNsInDisp(extrapOpts,stressInt.s);
end