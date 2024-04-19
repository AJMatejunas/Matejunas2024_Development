function [dispCrop,posCropX,posCropY] = func_cropDispFields(extrapOptsD,pos,disp)
% Crops bad edge data from the grid method or DIC data, used prior to used
% prior to extrapolating the edge data.
%
% Author: Lloyd Fletcher
% Image-Based Mechanics Group (IBMG), University of Southampton
% Date Created: 28/7/2020

    if strcmp(extrapOptsD.fieldOpt,'XY')
        % 'XY' option only crops in the direction pixels are lost (X for X
        % direction) - legacy option, good for symmetric data with no rotation
        dispCrop.x = func_cropFields_v4(disp.x,'X',...
            extrapOptsD.cropPx1st,extrapOptsD.cropPx1st);  
        posCropX.xGrid = func_cropFields_v4(pos.xGrid,'X',...
            extrapOptsD.cropPx1st,extrapOptsD.cropPx1st);
        posCropX.yGrid = func_cropFields_v4(pos.yGrid,'Y',...
            extrapOptsD.cropPx1st,extrapOptsD.cropPx1st); 

        dispCrop.y = func_cropFields_v4(disp.y,'Y',...
            extrapOptsD.cropPx1st,extrapOptsD.cropPx1st);
        posCropY = posCropX;
    elseif strcmp(extrapOptsD.fieldOpt,'both')
        % 'both' option crops in X and Y for both fields
        dispCrop.x = func_cropFields_v4(disp.x,'both',...
            extrapOptsD.cropPx1st,extrapOptsD.cropPx2nd);
        posCropX.xGrid = func_cropFields_v4(pos.xGrid,'both',...
            extrapOptsD.cropPx1st,extrapOptsD.cropPx2nd);
        posCropX.yGrid = func_cropFields_v4(pos.yGrid,'both',...
            extrapOptsD.cropPx1st,extrapOptsD.cropPx2nd);

        dispCrop.y = func_cropFields_v4(disp.y,'both',extrapOptsD.cropPx2nd,extrapOptsD.cropPx1st);
        posCropY.xGrid = func_cropFields_v4(pos.xGrid,'both',...
            extrapOptsD.cropPx2nd,extrapOptsD.cropPx1st);
        posCropY.yGrid = func_cropFields_v4(pos.yGrid,'both',...
            extrapOptsD.cropPx2nd,extrapOptsD.cropPx1st);
    end 
    
    posCropX.x = posCropX.xGrid(1,:);
    posCropX.y = flip(posCropX.yGrid(:,1)');
    posCropY.x = posCropY.xGrid(1,:);
    posCropY.y = flip(posCropY.yGrid(:,1)');
    
end

