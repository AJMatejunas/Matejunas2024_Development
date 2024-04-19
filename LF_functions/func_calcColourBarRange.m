function cRange = func_calcColourBarRange(calcType,plotVar,axisKernal,ampFactor)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 13/3/2017
%
% Calculates colourbar range based on the specified method

if nargin < 3
    if ~strcmp(calcType,'MaxQuantileOverCompOnly')
        % Best for in-plane CF with closing cracks
        ampFactor = 1.5;
        axisKernal = [0.02,0.98];

        % Through thickness CFRP with many cracks
        ampFactor = 1.5;
        axisKernal = [0.05,0.95];
    else
        ampFactor = 1.1;
        axisKernal = [0.02,0.98];    
    end
end

if strcmp(calcType,'Quantiles')
    % Get the quantile bounds
    QTs = quantile(plotVar(:),axisKernal);
    
    % Amplify the identified quantiles so we don't max the colorbar
    upperLim = round(ampFactor*QTs(2),1,'significant');
    lowerLim = round(ampFactor*QTs(1),1,'significant');   
    cRange = [lowerLim,upperLim];
    
elseif strcmp(calcType,'MaxQuantile')
    % Get the quantile bounds
    QTs = quantile(plotVar(:),axisKernal);
    maxQT = max(abs(QTs));
    
    % Amplify the identified quantiles so we don't max the colorbar
    upperLim = round(ampFactor*maxQT,2,'significant');
    lowerLim = round(-ampFactor*maxQT,2,'significant');
    cRange = [lowerLim,upperLim];
    
elseif strcmp(calcType,'MaxQuantileOverCompOnly')
    compRange = 1:ceil(size(plotVar,2)/2);
    tempVar = plotVar(:,compRange);
   % Get the quantile bounds
    QTs = quantile(tempVar(:),axisKernal);
    maxQT = max(abs(QTs));
    
    % Amplify the identified quantiles so we don't max the colorbar
    upperLim = round(ampFactor*maxQT,2,'significant');
    lowerLim = round(-ampFactor*maxQT,2,'significant');
    cRange = [lowerLim,upperLim];
       
elseif strcmp(calcType,'CentreMax')
    % Specify the centre portion of the specimen as ranges
    testRangeX = floor(size(plotVar,2)*0.25):ceil(size(plotVar,2)*0.75);
    testRangeY = floor(size(plotVar,1)*0.25):ceil(size(plotVar,1)*0.75);
    % Crop the plot variable based on the middle 50% range
    testVar = plotVar(testRangeY,testRangeX,:);
    % Get the min and max value in the range
    cRange = [min(testVar(:)),max(testVar(:))];
    clear testVar
else
    cRange = [0,0];
end 


end

