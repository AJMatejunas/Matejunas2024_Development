function [nX,nY] = func_getNormCoords(hax,Xin,Yin)
% Gets normalised axis co-ords for a given point on set of axes
xLims = xlim(hax);
yLims = ylim(hax);
xRange = max(xLims) - min(xLims);
yRange = max(yLims) - min(yLims);

nX = (Xin - min(xLims)) / xRange;
nY = (Yin - min(yLims)) / yRange;
end

