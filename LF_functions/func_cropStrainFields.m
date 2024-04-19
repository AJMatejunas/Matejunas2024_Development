function [strainCrop,posCropX,posCropY,posCropS] = func_cropStrainFields(...
    extrapOptsS,pos,strain)
% Crops bad edge data from the grid method or DIC data, used prior to
% extrapolating the edge data.
%
% Author: Lloyd Fletcher
% Image-Based Mechanics Group (IBMG), University of Southampton
% Date Created: 29/7/2020
% Date Edited: 31/7/2020

    if strcmp(extrapOptsS.fieldOpt,'XY')
        % 'XY' option only crops in the direction pixels are lost (X for X
        % direction) - legacy option, good for symmetric data with no rotation
        strainCrop.x = func_cropFields_v4(strain.x,'X',...
            extrapOptsS.extrapPx1st,extrapOptsS.extrapPx1st);  
        posCropX.xGrid = func_cropFields_v4(pos.xGrid,'X',...
            extrapOptsS.extrapPx1st,extrapOptsS.extrapPx1st);
        posCropX.yGrid = func_cropFields_v4(pos.yGrid,'Y',...
            extrapOptsS.extrapPx1st,extrapOptsS.extrapPx1st); 

        strainCrop.y = func_cropFields_v4(strain.y,'Y',...
            extrapOptsS.extrapPx1st,extrapOptsS.extrapPx1st);
        posCropY = posCropX;

    elseif strcmp(extrapOptsS.fieldOpt,'both')
        % 'both' option crops in X and Y for both fields
        strainCrop.x = func_cropFields_v4(strain.x,'both',...
            extrapOptsS.extrapPx1st,extrapOptsS.extrapPx2nd);
        posCropX.xGrid = func_cropFields_v4(pos.xGrid,'both',...
            extrapOptsS.extrapPx1st,extrapOptsS.extrapPx2nd);
        posCropX.yGrid = func_cropFields_v4(pos.yGrid,'both',...
            extrapOptsS.extrapPx1st,extrapOptsS.extrapPx2nd);

        strainCrop.y = func_cropFields_v4(strain.y,'both',...
            extrapOptsS.extrapPx2nd,extrapOptsS.extrapPx1st);
        posCropY.xGrid = func_cropFields_v4(pos.xGrid,'both',...
            extrapOptsS.extrapPx2nd,extrapOptsS.extrapPx1st);
        posCropY.yGrid = func_cropFields_v4(pos.yGrid,'both',...
            extrapOptsS.extrapPx2nd,extrapOptsS.extrapPx1st);     
    end
    
    % The shear strain is a 'cross' term so it is corrupted in both
    % directions, therefore primary crop in both directions regardless.
    strainCrop.s = func_cropFields_v4(strain.s,'both',...
        extrapOptsS.extrapPx1st,extrapOptsS.extrapPx1st);
    posCropS.xGrid = func_cropFields_v4(pos.xGrid,'both',...
        extrapOptsS.extrapPx1st,extrapOptsS.extrapPx1st);
    posCropS.yGrid = func_cropFields_v4(pos.yGrid,'both',...
        extrapOptsS.extrapPx1st,extrapOptsS.extrapPx1st);
    
    posCropX.x = posCropX.xGrid(1,:);
    posCropX.y = flip(posCropX.yGrid(:,1)');
    posCropY.x = posCropY.xGrid(1,:);
    posCropY.y = flip(posCropY.yGrid(:,1)');
    posCropS.x = posCropY.xGrid(1,:);
    posCropS.y = flip(posCropY.yGrid(:,1)');    
end

