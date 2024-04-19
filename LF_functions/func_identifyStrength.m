function [identStrSG,identStrLSG,identStrFS] = func_identifyStrength(...
    time,virtualGauge,fracture,stress,strain,strainRate,stiffness,method)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 27/3/2017
% Gets the fracture frame from user input and uses this to calculate the
% strength with the stress gauge and stress reconstructed from the strains

    %----------------------------------------------------------------------
    % Identify the Strength Using the Stress Gauge Approach 'Average'
    
    % Set the virtual gauge location to half a gauge width either side of
    % the selected point
    identStrSG.xRangePx = virtualGauge.xRange; 
    
    % Get the stress gauge values over the selected spatial range at the
    % specified fracture frame
    identStrSG.stressOverXRange = stress.xAvg(identStrSG.xRangePx,fracture.strengthFrame);
    identStrSG.stressOverXRange = identStrSG.stressOverXRange'; 
    % Take the average of the gauge locations to determine the strength
    identStrSG.strengthAvgOverXRange = mean(identStrSG.stressOverXRange);
    
    % Get the maximum tensile stress from the given stress gauge in the virtual gauge area 
    stressInGaugeArea = stress.xAvg(identStrSG.xRangePx,:);
    [identStrSG.maxTensStressOverXRange,ind] = max(stressInGaugeArea(:));
    [maxTensLocXInGaugeArea,maxTensFrameInGaugeArea] = ...
        ind2sub(size(stressInGaugeArea),ind);

    % Get the strength from the user defined fracture location and frame
    identStrSG.strengthAtFractLocX = stress.xAvg(fracture.locX,fracture.strengthFrame);
    
    % Get the maximum tensile stress at the fracture location over all time
    [identStrSG.maxTensStressAtFracLocX,identStrSG.maxTensFrameAtFracLocX] =...
        max(stress.xAvg(fracture.locX,:));

    %----------------------------------------------------------------------
    % Identify the Strength using the Linear Stress Gauge
    identStrLSG.xRangePx = virtualGauge.xRange;
    identStrLSG.yRangePx = virtualGauge.yRange;

    % Get the stress gauge values over the selected spatial range at the
    % specified fracture frame
    identStrLSG.stressInGaugeAreaAtFracFrame = stress.xLinearGauge(...
        virtualGauge.yRange,virtualGauge.xRange,fracture.strengthFrame);
    % Take the average of the gauge locations to determine the strength
    identStrLSG.strengthAvgOverGaugeArea = mean(mean(identStrLSG.stressInGaugeAreaAtFracFrame));
    
    % Get the maximum tensile stress from the given stress gauge in the virtual gauge area 
    stressInGaugeArea = stress.xLinearGauge(virtualGauge.yRange,virtualGauge.xRange,:);
    [identStrLSG.maxTensStressInGaugeArea,ind] = max(stressInGaugeArea(:));
    [maxTensLocYInGaugeArea,maxTensLocXInGaugeArea,...
        maxTensFrameInGaugeArea] = ind2sub(size(stressInGaugeArea),ind);
    
    % Get the average stress over the gauge line at frac loc X
    identStrLSG.strengthAvgGaugeYFracLineX = mean(stress.xLinearGauge(...
        virtualGauge.yRange,fracture.locX,fracture.strengthFrame));
    
    % Get the average stress over the gauge line at frac loc X
    identStrLSG.maxTensStressGaugeYFracLineX = max(stress.xLinearGauge(...
        virtualGauge.yRange,fracture.locX,fracture.strengthFrame));

    % Get the strength from the user defined fracture location and frame
    identStrLSG.strengthAtFractLocXY = stress.xLinearGauge(...
        fracture.locY,fracture.locX,fracture.strengthFrame);
    
    % Get the maximum tensile stress at the fracture location over all time
    [identStrLSG.maxTensStressAtFracLocXY,identStrLSG.maxTensFrameAtFracLocXY] =...
        max(stress.xLinearGauge(fracture.locY,fracture.locX,:));
   
    %----------------------------------------------------------------------
    % Identify the Strength Using the Reconstructed Stress Field 'Local'
    if stiffness.method == 1
        stressQRecon = stiffness.Ex*strain.x;
    else
        stressQRecon = stiffness.Ex/(1-stiffness.nuxy^2)*(strain.x+stiffness.nuxy*strain.y);
    end 
    % Use the virtual gauge area to calculate the failure strain
    identStrFS.xRangePx = virtualGauge.xRange;
    identStrFS.yRangePx = virtualGauge.yRange;
    % Calculate the failure strain over the gauge area
    identStrFS.avgStrainInGaugeArea = squeeze(mean(mean(strain.x(virtualGauge.yRange,virtualGauge.xRange,fracture.strengthFrame))));
    identStrFS.strengthFromGaugeAvg = squeeze(mean(mean(stressQRecon(virtualGauge.yRange,virtualGauge.xRange,fracture.strengthFrame))));
    % Calculate a pseudo stress-gauge equivalent from the strains
    identStrFS.avgAxialStrain = squeeze(mean(strain.x(:,fracture.locX,fracture.strengthFrame)));
    identStrFS.stressFromAvgAxialStrain = squeeze(mean(stressQRecon(:,fracture.locX,fracture.strengthFrame)));    
    
end

