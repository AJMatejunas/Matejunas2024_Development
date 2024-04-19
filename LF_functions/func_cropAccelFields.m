function [accelCrop,posCropX,posCropY] = func_cropAccelFields(...
    extrapOptsA,pos,accel)
% Crops bad edge data from the grid method or DIC data, used prior to 
% extrapolating the edge data.
%
% Author: Lloyd Fletcher
% Image-Based Mechanics Group (IBMG), University of Southampton
% Date Created: 29/7/2020
% Date Edited: 31/7/2020

    if strcmp(extrapOptsA.fieldOpt,'XY')
        % 'XY' option only crops in the direction pixels are lost (X for X
        % direction) - legacy option, good for symmetric data with no rotation
        accelCrop.x = func_cropFields_v4(accel.x,'X',...
            extrapOptsA.extrapPx1st,extrapOptsA.extrapPx1st);  
        posCropX.xGrid = func_cropFields_v4(pos.xGrid,'X',...
            extrapOptsA.extrapPx1st,extrapOptsA.extrapPx1st);
        posCropX.yGrid = func_cropFields_v4(pos.yGrid,'Y',...
            extrapOptsA.extrapPx1st,extrapOptsA.extrapPx1st); 

        accelCrop.y = func_cropFields_v4(accel.y,'Y',...
            extrapOptsA.extrapPx1st,extrapOptsA.extrapPx1st);
        posCropY = posCropX;

    elseif strcmp(extrapOptsA.fieldOpt,'both')
        % 'both' option crops in X and Y for both fields
        accelCrop.x = func_cropFields_v4(accel.x,'both',...
            extrapOptsA.extrapPx1st,extrapOptsA.extrapPx2nd);
        posCropX.xGrid = func_cropFields_v4(pos.xGrid,'both',...
            extrapOptsA.extrapPx1st,extrapOptsA.extrapPx2nd);
        posCropX.yGrid = func_cropFields_v4(pos.yGrid,'both',...
            extrapOptsA.extrapPx1st,extrapOptsA.extrapPx2nd);

        accelCrop.y = func_cropFields_v4(accel.y,'both',...
            extrapOptsA.extrapPx2nd,extrapOptsA.extrapPx1st);
        posCropY.xGrid = func_cropFields_v4(pos.xGrid,'both',...
            extrapOptsA.extrapPx2nd,extrapOptsA.extrapPx1st);
        posCropY.yGrid = func_cropFields_v4(pos.yGrid,'both',...
            extrapOptsA.extrapPx2nd,extrapOptsA.extrapPx1st);     
    end
   
    posCropX.x = posCropX.xGrid(1,:);
    posCropX.y = flip(posCropX.yGrid(:,1)');
    posCropY.x = posCropY.xGrid(1,:);
    posCropY.y = flip(posCropY.yGrid(:,1)');  
end

