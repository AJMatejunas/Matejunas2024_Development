function [identStiffness,linearFitCoeffs,fitFrameRange] = ...
    func_identStiffLinFitStressStrainCurve(identOpts,stressAvg,strainAvg)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 9/2/2018
% Date Modified: 20/9/2018 
%
% Linearly fits the stressAvg vs strainAvg to obtain the stiffness

% Turn off warnings:
% Stops many warnings being thrown when the linear fit is bad.
warning('off','all')

%--------------------------------------------------------------------------
% Specify the range over which the fit is performed

% Deafult starting index given by the specified range
startInd = identOpts.fitDefaultRange(1);

if strcmp(identOpts.limParam,'strain')
    thresholdParam = strainAvg;
else
    thresholdParam = stressAvg;
end

% Option 1: fit the curve only over the compressive loading
if strcmp(identOpts.fitRangeOpt,'componly')
    for xx = 1:size(thresholdParam,1)
        [~,minInd] = min(thresholdParam(xx,:));
        % Allows the start point of the fit to be thresholded above a
        % certain value
        startInd = identOpts.fitDefaultRange(1);
        if identOpts.minThreshold ~= 0
            for tt = 1:minInd
                if abs(thresholdParam(xx,tt)) > abs(identOpts.minThreshold) 
                    startInd = tt;
                    break;
                end
            end
        end
        % Store the fit range in a cell vector because the range will vary
        % based on position
        fitFrameRange{xx} = startInd:minInd;
    end
    
% Option 2: fit the curve only over the tensile loading
elseif strcmp(identOpts.fitRangeOpt,'tensonly')
    for xx = 1:size(thresholdParam,1)
        [~,maxInd] = max(thresholdParam(xx,:));
        % Allows the start point of the fit to be thresholded above a
        % certain value
        startInd = identOpts.fitDefaultRange(1);
        if identOpts.minThreshold ~= 0
            for tt = 1:maxInd
                if abs(thresholdParam(xx,tt)) > abs(identOpts.minThreshold) 
                    startInd = tt;
                    break;
                end
            end
        end
        fitFrameRange{xx} = startInd:maxInd;
    end
    
% Option 3: fit the curve over the compressive or tensile portion of the
% test whichever has the highest maximum stress
elseif strcmp(identOpts.fitRangeOpt,'maxabs')
    for xx = 1:size(thresholdParam,1)
        [minS,minInd] = min(thresholdParam(xx,:));
        [maxS,maxInd] = max(thresholdParam(xx,:));

        if abs(minS) > abs(maxS)
             fitFrameRange{xx} = startInd:minInd;
        else
             fitFrameRange{xx} = startInd:maxInd;
        end
    end
    
% Option 4: fit the curve up to the specified threshold value
elseif strcmp(identOpts.fitRangeOpt,'maxThreshold')
    for xx = 1:size(thresholdParam,1)
        % If there is a minimum threshold find the starting frame for the fit
        % for each slice
        startInd(xx) = identOpts.fitDefaultRange(1);
        if identOpts.minThreshold ~= 0
            for tt = 1:size(thresholdParam,2)
                if abs(thresholdParam(xx,tt)) > abs(identOpts.minThreshold) 
                    startInd(xx) = tt;
                    break;
                end
            end
        end
    end
    
    for xx = 1:size(thresholdParam,1)
        for tt = 1:size(thresholdParam,2)
             endInd(xx) = size(thresholdParam,2);
             if abs(thresholdParam(xx,tt)) > abs(identOpts.maxThreshold) 
                endInd(xx) = tt-1;
                break;
            end  
        end
        fitFrameRange{xx} = startInd(xx):endInd(xx);
    end
     
% Default: fit the curve over the specified default time range  
else
    for xx = 1:size(thresholdParam,1)
        % Allows the start point of the fit to be thresholded above a
        % certain value
        startInd = identOpts.fitDefaultRange(1);
        if identOpts.minThreshold ~= 0
            for tt = 1:size(thresholdParam,2)
                if abs(thresholdParam(xx,tt)) > abs(identOpts.minThreshold) 
                    startInd = tt;
                    break;
                end
            end
        end
        fitFrameRange{xx} = startInd:identOpts.fitDefaultRange(end);
    end
end

%--------------------------------------------------------------------------
% Specify the fit options
if identOpts.robustFit
    options = fitoptions('Method','LinearLeastSquares','Robust','On');
else
    options = fitoptions('Method','LinearLeastSquares','Robust','Off');
end

%--------------------------------------------------------------------------
% Loop over each transverse slice and fit the stress strain curves

% Pre-alloc Qxx vars for speed
QFit = cell(1,size(stressAvg,1));
QGOF = cell(1,size(stressAvg,1));
identStiffness = zeros(1,size(stressAvg,1));

for xx = 1:size(stressAvg,1)
    % Make sure there are enough points to linearly fit
    if (length(fitFrameRange{xx}) >= 2)
        if (sum(isnan(strainAvg(xx,fitFrameRange{xx}))) == 0) && (sum(isnan(stressAvg(xx,fitFrameRange{xx}))) == 0)

           [QFit{xx},QGOF{xx}] = fit(strainAvg(xx,fitFrameRange{xx})',...
               stressAvg(xx,fitFrameRange{xx})','poly1',options);

           identStiffness(xx) = QFit{xx}.p1;
           linearFitCoeffs{xx} = [QFit{xx}.p1,QFit{xx}.p2]; 
        else
            fprintf('WARNING: NaNs in stress-strain curve at pixel %i.\n',xx)
            identStiffness(xx) = NaN;
            linearFitCoeffs{xx} = [NaN,NaN]; 
        end
    else
       % If there are not enough points return 0 stiffness
       % fprintf('WARNING: not enough points for a linear fit of stress-strain curve.\n')
       identStiffness(xx) = 0;
       linearFitCoeffs{xx} = [0,0];
    end
end

% Turn on warnings:
warning('on','all')
    
end

