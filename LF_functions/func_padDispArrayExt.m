function [dispTemp] = func_padDispArrayExt(pos,dispCrop,grid,padRadius,padEdgeMethod,cropPx,fitWindow,extrapOpts) 
% This function pads out ux and uy displacement fields to recover lost
% information at the edges. The data is padded out by 'padRadius' at the
% edges. This padding is performed in three ways based on the option
% specified in 'extrapOpts.padEdgeMethod' ('quadratic', 'linear
% regression', or 'linear'(two point extrapolation - interp1). The padded 
% displacement fields are stored in 'dispTemp'

% PhotoDyn Group
% Author: Jared Van Blitterswyk
% Modified By: Lloyd Fletcher
% Last modified: 2/5/2018

[sy,sx,st] = size(dispCrop.x);
[syY,sxY,~] = size(dispCrop.y);

if isfield(grid,'asymmPitch')
    if grid.asymmPitch
        % define cropped and padded x and y coordinates
        pos.xCrop = pos.x(cropPx+1:end-cropPx);
        pos.xPad = min(pos.xCrop)-padRadius*grid.mPerPxX:grid.mPerPxX:max(pos.xCrop)+padRadius*grid.mPerPxX;

        pos.yCrop = pos.y(cropPx+1:end-cropPx);
        pos.yPad = min(pos.yCrop)-padRadius*grid.mPerPxY:grid.mPerPxY:max(pos.yCrop)+padRadius*grid.mPerPxY;
    else
        % define cropped and padded x and y coordinates
        pos.xCrop = pos.x(cropPx+1:end-cropPx);
        pos.xPad = min(pos.xCrop)-padRadius*grid.mPerPx:grid.mPerPx:max(pos.xCrop)+padRadius*grid.mPerPx;

        pos.yCrop = pos.y(cropPx+1:end-cropPx);
        pos.yPad = min(pos.yCrop)-padRadius*grid.mPerPx:grid.mPerPx:max(pos.yCrop)+padRadius*grid.mPerPx;
    end
else
    % define cropped and padded x and y coordinates
    pos.xCrop = pos.x(cropPx+1:end-cropPx);
    pos.xPad = min(pos.xCrop)-padRadius*grid.mPerPx:grid.mPerPx:max(pos.xCrop)+padRadius*grid.mPerPx;

    pos.yCrop = pos.y(cropPx+1:end-cropPx);
    pos.yPad = min(pos.yCrop)-padRadius*grid.mPerPx:grid.mPerPx:max(pos.yCrop)+padRadius*grid.mPerPx;
end

if extrapOpts.imageDef == true
    dispTemp.x = zeros(sy+2*padRadius,sx+2*padRadius,st);
    dispTemp.x(padRadius+1:end-padRadius,padRadius+1:end-padRadius,:) = dispCrop.x;

    dispTemp.y = zeros(syY+2*padRadius,sxY+2*padRadius,st);
    dispTemp.y(padRadius+1:end-padRadius,padRadius+1:end-padRadius,:) = dispCrop.y;    
else
    dispTemp.x = zeros(sy,sx+2*padRadius,st);
    dispTemp.x(:,padRadius+1:end-padRadius,:) = dispCrop.x;

    dispTemp.y = zeros(syY+2*padRadius,sxY,st);
    dispTemp.y(padRadius+1:end-padRadius,:,:) = dispCrop.y;
end

%% 1D interpolation
% create vectors of data to fit
xToFit1 = pos.xCrop(1+grid.pxPerPeriod:1+fitWindow+grid.pxPerPeriod);
xToFit2 = pos.xCrop(end-grid.pxPerPeriod-fitWindow:end-grid.pxPerPeriod);

xToExtrap1 = pos.xPad(1:padRadius);
xToExtrap2 = pos.xPad(end-padRadius:end);

for t = 1:st
for p = 1:sy   
    % create vectors of displacements to fit
    % fit smoothed map from edges inward by fitWindow plus one pitch
    dispXToFit1 = dispCrop.x(p,1+grid.pxPerPeriod:1+fitWindow+grid.pxPerPeriod,t);
    dispXToFit2 = dispCrop.x(p,end-fitWindow-grid.pxPerPeriod:end-grid.pxPerPeriod,t);
    
    if strcmp(padEdgeMethod,'quadratic')
        % perform fit to cropped field     
        x1 = [ones(length(xToFit1),1) xToFit1' xToFit1'.^2];
        x2 = [ones(length(xToFit2),1) xToFit2' xToFit2'.^2];
        r1 = x1\dispXToFit1';
        r2 = x2\dispXToFit2';
        
        % evaluate function at padded coordinates
        dispTemp.x(p,1:padRadius,t) = r1(3)*xToExtrap1.^2 + r1(2)*xToExtrap1 + r1(1);
        dispTemp.x(p,end-padRadius:end,t) = r2(3)*xToExtrap2.^2 + r2(2)*xToExtrap2 + r2(1);
        
    elseif strcmp(padEdgeMethod,'regression')
        % perform fit to cropped field
        x1 = [ones(length(xToFit1),1) xToFit1'];
        x2 = [ones(length(xToFit2),1) xToFit2'];
        r1 = x1\dispXToFit1';
        r2 = x2\dispXToFit2';
        
        diff1 = dispTemp.x(p,padRadius+1,t) - (r1(2)*pos.xCrop(1) + r1(1));
        diff2 = dispTemp.x(p,end-padRadius,t) - (r2(2)*pos.xCrop(end) + r2(1));
        
        r1(1) = r1(1) + diff1;
        r2(1) = r2(1) + diff2;
        
        % evaluate function at padded coordinates
        dispTemp.x(p,1:padRadius,t) = r1(2)*xToExtrap1 + r1(1);
        dispTemp.x(p,end-padRadius:end,t) = r2(2)*xToExtrap2 + r2(1);
        
    else % default MATLAB extrapolation
        % displacement fields to perform padding
        temp_extrapVar(p,:,t) = interp1(pos.xCrop,dispCrop.x(p,:,t),pos.xPad,padEdgeMethod,'extrap');

        % Add the extrapolated data back onto the edges
        dispTemp.x(p,1:padRadius,t) = temp_extrapVar(p,1:padRadius,t);
        dispTemp.x(p,end-padRadius+1:end,t) = temp_extrapVar(p,end-padRadius+1:end,t);
    end
    if extrapOpts.imageDef == true % data cropped on all borders, need to pad out ux as well
        temp_extrapVar2(p,:,t) = interp1(pos.xCrop,dispCrop.y(p,:,t),pos.xPad,'linear','extrap');
        % Add the extrapolated data back onto the edges
        dispTemp.y(p,1:padRadius,t) = temp_extrapVar2(p,1:padRadius,t);
        dispTemp.y(p,end-padRadius+1:end,t) = temp_extrapVar2(p,end-padRadius+1:end,t);
    end
end

end
 clearvars temp_extrapVar
 if extrapOpts.imageDef == true
     clearvars temp_extrapVar2
 end
 
yToFit1 = pos.yCrop(1+grid.pxPerPeriod:1+fitWindow+grid.pxPerPeriod);
yToFit2 = pos.yCrop(end-grid.pxPerPeriod-fitWindow:end-grid.pxPerPeriod);

yToExtrap1 = pos.yPad(1:padRadius);
yToExtrap2 = pos.yPad(end-padRadius:end);
 
for t = 1:st
for p = 1:sxY
    
    % perform fit to cropped field
    % fit smoothed map from edges inward by fitWindow plus one pitch
    dispYToFit1 = dispCrop.y(1+grid.pxPerPeriod:1+fitWindow+grid.pxPerPeriod,p,t);
    dispYToFit2 = dispCrop.y(end-fitWindow-grid.pxPerPeriod:end-grid.pxPerPeriod,p,t);
    
    if strcmp(padEdgeMethod,'quadratic')
        % perform fit to cropped field
        
        y1 = [ones(length(yToFit1),1) yToFit1' yToFit1'.^2];
        y2 = [ones(length(yToFit2),1) yToFit2' yToFit2'.^2];
        r1 = y1\dispYToFit1;
        r2 = y2\dispYToFit2;
        
        % evaluate function at padded coordinates
        dispTemp.y(1:padRadius,p,t) = r1(3)*yToExtrap1.^2 + r1(2)*yToExtrap1 + r1(1);
        dispTemp.y(end-padRadius:end,p,t) = r2(3)*yToExtrap2.^2 + r2(2)*yToExtrap2 + r2(1);
        
    elseif strcmp(padEdgeMethod,'regression')
        % perform fit to cropped field       
        y1 = [ones(length(yToFit1),1) yToFit1'];
        y2 = [ones(length(yToFit2),1) yToFit2'];
        r1 = y1\dispYToFit1;
        r2 = y2\dispYToFit2;
        
        diff1 = dispTemp.y(padRadius+1,p,t) - (r1(2)*pos.yCrop(1) + r1(1));
        diff2 = dispTemp.y(end-padRadius,p,t) - (r2(2)*pos.yCrop(end) + r2(1));
        
        r1(1) = r1(1) + diff1;
        r2(1) = r2(1) + diff2;
        
        % evaluate function at padded coordinates
        dispTemp.y(1:padRadius,p,t) = r1(2)*yToExtrap1 + r1(1);%dispTemp.x(p,padRadius+1,t);
        dispTemp.y(end-padRadius:end,p,t) = r2(2)*yToExtrap2 + r2(1);%dispTemp.x(p,end-padRadius,t);
        
    else  % default MATLAB extrapolation
        % displacement fields to perform padding
        temp_extrapVar(:,p,t) = interp1(pos.yCrop,dispCrop.y(:,p,t),pos.yPad,padEdgeMethod,'extrap');
        % Add the extrapolated data back onto the edges
        dispTemp.y(1:padRadius,p,t) = temp_extrapVar(1:padRadius,p,t);
        dispTemp.y(end-padRadius+1:end,p,t) = temp_extrapVar(end-padRadius+1:end,p,t);
    end
    if extrapOpts.imageDef == true % data cropped on all borders, need to pad out ux as well
        temp_extrapVar2(:,p,t) = interp1(pos.yCrop,dispCrop.x(:,p,t),pos.yPad,'linear','extrap');
        % Add the extrapolated data back onto the edges
        dispTemp.x(1:padRadius,p,t) = temp_extrapVar2(1:padRadius,p,t);
        dispTemp.x(end-padRadius+1:end,p,t) = temp_extrapVar2(end-padRadius+1:end,p,t);
    end
end
end