function [mapZ,minErr] = func_contourImage3VarsSubPlot(plotParams,xVar,yVar,zVar)
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

% Upsamples the data by linear interpolation
zInterp = scatteredInterpolant(yVar,xVar,zVar);
xUS = linspace(min(X),max(X),plotParams.imageUpSample);
yUS = linspace(min(Y),max(Y),plotParams.imageUpSample);
[xUSGrid,yUSGrid] = meshgrid(xUS,yUS);
zUS = zInterp(yUSGrid,xUSGrid);

% Find the minimum z value and its index
[minVal,ind] = min(abs(mapZ(:)));
[minSK,minTK] = ind2sub(size(mapZ),ind);
minErr = [yVals(minSK),xVals(minTK),minVal,minSK,minTK];

% Calculate the colourbar range
errMin = min(mapZ(:));
errMax = max(mapZ(:));

if abs(errMax) > abs(errMin)
    cBounds = [errMin,quantile(mapZ(:),0.9)];
else
    cBounds = [quantile(mapZ(:),0.1),errMax];
end

% Plot the contour map
imagesc(yUS,xUS,zUS)
colorbar
colormap(jet)
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
set(gca,'YDir','normal')

end

