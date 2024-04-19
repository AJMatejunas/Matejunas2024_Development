function func_plotSweepData(plotParams,data,dataHeaders,dataUnits,plotConfig,labelStrs)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 10/4/2017
%
%   Plots parametric sweep data 
%   Expects the following inputs:
%   data: matrix of sweep variables arranged column wise
%   dataHeaders: cell array of strings giving headers for columns of data
%   dataUnits: cell array of strings as units for each column of data
%   plotConfig: length 5 vector giving indices for the following
%        plotConfig(1) = column of x axis variable
%        plotConfig(2) = column of y axis variable
%        plotConfig(3) = column of variable to separate plot lines
%        plotConfig(4) = column of fixed variable
%        plotConfig(5) = value of fixed variable      

% Assign variables to required indices
xInd = plotConfig(1);
yInd = plotConfig(2);
lInd = plotConfig(3);
fInd = plotConfig(4);
fVal = plotConfig(5);

% Check if the label strings have been specified otherwise use defaults
if nargin < 6
    labelStrs.x = [dataHeaders{xInd},' ',dataUnits{xInd}];
    labelStrs.y = [dataHeaders{yInd},' ',dataUnits{yInd}];
    if fInd > 0
        labelStrs.t = [dataHeaders{fInd},'= ',num2str(fVal)];
    end
end

% Make sure the required fields are in the plot params struct
if ~isfield(plotParams,'legLoc')
    plotParams.legLoc = 'best';
end
if ~isfield(plotParams,'titleOn')
    plotParams.titleOn = true;
end
if ~isfield(plotParams,'formatType')
    plotParams.formatType = 'diagnostic';
end
if ~isfield(plotParams,'interpret')
    plotParams.interpret = 'tex';
end


% Create the plot props structure for the correct formatting
plotProps = func_initPlotPropsStruct(plotParams.formatType);

% Sort the data to get only the data for the fixed param
[numRows,numCols] = size(data);
% if the fixed index is zero we create a dummy column and use all the data
thirdVarFlag = true;
if fInd == 0
    thirdVarFlag = false;
    fixedData = [data,zeros(numRows,1)];
    find = numCols+1;
    fval = 0;
else
    sortedData = sortrows(data,fInd);
    mask = sortedData(:,fInd) == fVal;
    fixedData = sortedData(mask,:);
end

% Sort the data by the line parameter
fixedData = sortrows(fixedData,lInd);
% find the unique vals and the number of them
uniqueVals = unique(fixedData(:,lInd));
uniqueCount = [uniqueVals,histc(fixedData(:,lInd),uniqueVals)];

startI = 1;
endI = uniqueCount(1,2);
for i = 1:length(uniqueVals)
   plotX{i} = fixedData(startI:endI,xInd);
   plotY{i} = fixedData(startI:endI,yInd);
   [minY,indMinY] = min(fixedData(startI:endI,yInd));
   minX = fixedData((startI+indMinY-1),xInd);
   minPt{i} = [minX,minY];
   legendStrs{i} = [dataHeaders{lInd},'=',num2str(uniqueVals(i)),dataUnits{lInd}];
   if i ~= length(uniqueVals) 
       startI = endI+1;
       endI = endI+uniqueCount(i,2);
   end
end

hold on
for i = 1:length(plotX)
   plot(plotX{i},plotY{i},plotProps.lineStyleVec{i},'linewidth',plotProps.lw,'markersize',plotProps.ms)     
end
hold off

if thirdVarFlag && plotParams.titleOn
    title(labelStrs.t,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'interpreter',plotParams.interpret)
end
xlabel(labelStrs.x,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'interpreter',plotParams.interpret)
ylabel(labelStrs.y,'fontsize',plotProps.hfs,'fontname',plotProps.ft,'interpreter',plotParams.interpret)
set(gca,'fontsize', plotProps.fs,'fontname',plotProps.ft)
set(gca,'XMinorTick','on','YMinorTick','on')
legend(legendStrs,'Location',plotParams.legLoc)
box on
grid on

end

