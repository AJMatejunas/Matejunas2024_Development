function [hf,mapZ,minErr,err00] = func_contour3Vars(plotParams,xVar,yVar,zVar)
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
minErr = [yVals(minSK),xVals(minTK),minVal];

% Create and size the figure
plotProps.sizePerFigXcm = plotProps.sizePerFigXcm*2;
plotProps.sizePerFigYcm = plotProps.sizePerFigYcm*2;
hf = func_createFigure(plotProps);

% Plot the contour map
contourf(X,Y,Z)
colorbar
colormap(jet)
if isfield(plotParams,'cAxis')
    caxis(plotParams.cAxis)
end

% Show the minimum value and it's location
if plotParams.plotMin
hold on
    plot(minErr(2),minErr(1),'+r','linewidth',2,'markersize',15)
hold off
end

title(plotParams.tStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'interpreter',plotParams.interpret)
xlabel(plotParams.xStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'interpreter',plotParams.interpret)
ylabel(plotParams.yStr,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'interpreter',plotParams.interpret)

end

