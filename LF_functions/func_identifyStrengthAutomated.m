function [identStrSG,identStrLSG,fracture] = func_identifyStrengthAutomated(...
    virtualGauge,fracture,stress)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 30/8/2017
% Identifies strength automatically taking the maximum stress gauge value

    % Specify the range over which the strength can be identified
    if isfield(fracture,'strengthFrameMax')
        stressFRange = 1:fracture.strengthFrameMax;
    else
        stressFRange = 1:size(stress.xAvg,2); 
    end

    %----------------------------------------------------------------------
    % Identify the Strength Using the Stress Gauge Approach 'Average'
    
    % Identify the 'max' strength over the virtual gauge 'X' area
    %{
    tempStress = stress.xAvg(virtualGauge.xRange,stressFRange );
    [identStrSG.rangeStrMax,tempInd] = max(tempStress(:));
    [identStrSG.rangeLocX,identStrSG.rangeStrFrame] = ind2sub(size(tempStress),tempInd);
    %}
    
    % Identify the 'point' strength at the exact specified fracture location
    tempStress = squeeze(stress.xAvg(fracture.locX,stressFRange));
    identStrSG.pointX = fracture.locX;
    [identStrSG.pointStrMax,identStrSG.pointStrFrame] = max(tempStress);
    
    % Identify the 'mean' strength over the virtual gauge 'X' area
    tempStress = mean(stress.xAvg(virtualGauge.xRange,stressFRange ));
    [identStrSG.virtGAvgStrMax,identStrSG.virtGAvgStrFrame] = max(tempStress);
    
    %----------------------------------------------------------------------
    % Identify the Strength using the Linear Stress Gauge
    
    % Identify the 'point' strength over the virtual gauge 'X' and 'Y' area
    %{
    tempStress = stress.xLinearGauge(virtualGauge.yRange,virtualGauge.xRange,stressFRange);
    [identStrLSG.rangeStrMax,tempInd] = max(tempStress(:));
    [identStrLSG.rangeLocY,identStrLSG.rangeLocX,identStrLSG.rangeStrFrame] = ...
        ind2sub(size(tempStress),tempInd);
    %}
    
    % Identify the 'point' strength at the exact specified fracture location
    tempStress = stress.xLinearGauge(fracture.locY,fracture.locX,stressFRange);
    [identStrLSG.pointStrMax,identStrLSG.pointStrFrame] = max(tempStress);
    
    % Identify the 'mean' strength over the virtual gauge 'X' and 'Y' area
    tempStress = mean(mean(stress.xLinearGauge(virtualGauge.yRange,virtualGauge.xRange,stressFRange)));
    [identStrLSG.virtGAvgStrMax,identStrLSG.virtGAvgStrFrame] = max(tempStress);
    
    %----------------------------------------------------------------------
    fracture.strengthFrame = identStrLSG.virtGAvgStrFrame;
    
end

