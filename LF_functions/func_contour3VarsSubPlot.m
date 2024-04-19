function [mapZ,minErr] = func_contour3VarsSubPlot(plotParams,xVar,yVar,zVar)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 24/1/2018
% Creates a contour plot from 3 input vectors

plotProps = func_initPlotPropsStruct(plotParams.formatType);

% Sort the z var into the 2D contour map
xVals = unique(xVar);
yVals = unique(yVar);
for i = 1:length(zVar)
    for x = 1:length(xVals)
        for y = 1:length(yVals)
            if xVar(i) == xVals(x)  && yVar(i) == yVals(y)
                mapZ(y,x) = zVar(i);
            end
        end
    end
end
Y = yVals;
X = xVals;
Z = mapZ;

% Find the minimum z value and its index
[minVal,ind] = min(abs(mapZ(:)));
[minSK,minTK] = ind2sub(size(mapZ),ind);
minErr = [yVals(minSK),xVals(minTK),minVal,minSK,minTK];

% Calculate the contours
xTemp = linspace(0,3.5,plotParams.numContours);
yTemp = 1./exp(xTemp);
errMin = min(mapZ(:));
errMax = max(mapZ(:));
contsEven = linspace(errMin,errMax,plotParams.numContours);

if abs(errMin) > abs(errMax)
    conts = yTemp.*contsEven;
    cBounds = [quantile(conts,0.25),errMax];
else
    conts = yTemp.*flip(contsEven);
    cBounds = [errMin,quantile(conts,0.75)];
end

% Plot the contour map
contourf(X,Y,Z,conts)
colorbar
colormap(jet)
caxis(cBounds)
if isfield(plotParams,'cAxis')
    caxis(plotParams.cAxis)
else
    caxis(cBounds)
end
%set(gca,'ColorScale','log')

% Show the minimum value and it's location
if plotParams.plotMin
hold on
    plot(minErr(2),minErr(1),'+r','linewidth',2,'markersize',15)
hold off
end

title(plotParams.tStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'interpreter',plotProps.interpText)
xlabel(plotParams.xStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'interpreter',plotProps.interpText)
ylabel(plotParams.yStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'interpreter',plotProps.interpText)
set(gca,'fontsize', plotProps.fs,'fontname',plotProps.ft)

end

