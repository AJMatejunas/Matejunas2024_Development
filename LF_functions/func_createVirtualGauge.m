function virtualGauge = func_createVirtualGauge(virtualGauge,fracture,strain)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 20/6/2017
%
% Defines the required bounds for the virtual gauge used for local strength
% identification. Ensures values out of the image are not used.


    virtualGauge.xRad = floor(virtualGauge.Xpx/2);
    
    virtualGauge.yRad = floor(virtualGauge.Ypx/2);
    virtualGauge.xRange = (fracture.locX-virtualGauge.xRad):(fracture.locX+virtualGauge.xRad);
    virtualGauge.yRange = (fracture.locY-virtualGauge.yRad):(fracture.locY+virtualGauge.yRad);
    
    % Check that our virtual strain gauge is within the bounds of the image
    % if not shift it
    
    % Check the X Range
    if virtualGauge.xRange(1) < 1
        pxCorrect = abs(1-virtualGauge.xRange(1));
        virtualGauge.xRange = virtualGauge.xRange + pxCorrect;
    end
    if virtualGauge.xRange(end) >  size(strain.x,2)
        pxCorrect = abs(virtualGauge.xRange(end)-size(strain.x,2));
        virtualGauge.xRange = virtualGauge.xRange - pxCorrect;    
    end
    % Check the Y Range
    if virtualGauge.yRange(1) < 1
        pxCorrect = abs(1-virtualGauge.yRange(1));
        virtualGauge.yRange = virtualGauge.yRange + pxCorrect;
    end
    if virtualGauge.yRange(end) >  size(strain.x,1)
        pxCorrect = abs(virtualGauge.yRange(end)-size(strain.x,1));
        virtualGauge.yRange = virtualGauge.yRange - pxCorrect;    
    end

end

