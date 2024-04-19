function [stress,strain] = func_calcAvgStressStrainInGaugeArea(globalOpts,strOpts,stress,strain)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 9/2/2018
%
% Calculates the mean stress/strain in the virtual gauge area for strength 
% strength identification purposes.Depending on the material model being 
% used the axial strain calculation will change.

virtualGauge = strOpts.virtualGauge;
if strcmp('isotropic',globalOpts.matModel)
    stress.virtGaugeAvg.x = squeeze(nanmean(nanmean(stress.xLinearGauge(virtualGauge.yRange,virtualGauge.xRange,:))));
    strain.virtGaugeAvg.x = squeeze(nanmean(nanmean(strain.x(virtualGauge.yRange,virtualGauge.xRange,:))));
    strain.virtGaugeAvg.y = squeeze(nanmean(nanmean(strain.y(virtualGauge.yRange,virtualGauge.xRange,:))));
    nuxy = strOpts.stiffnessQxy/strOpts.stiffnessQxx;
    strain.virtGaugeAvg.xny = strain.virtGaugeAvg.x + nuxy*strain.virtGaugeAvg.y;
    
    stress.QRecon.virtGaugeAvg.x = squeeze(nanmean(nanmean(stress.QRecon.x(virtualGauge.yRange,virtualGauge.xRange,:))));
  
elseif strcmp('orthotropicReduced',globalOpts.matModel)
    stress.virtGaugeAvg.x = squeeze(nanmean(nanmean(stress.xLinearGauge(virtualGauge.yRange,virtualGauge.xRange,:))));
    strain.virtGaugeAvg.x = squeeze(nanmean(nanmean(strain.x(virtualGauge.yRange,virtualGauge.xRange,:))));
     
    stress.QRecon.virtGaugeAvg.x = squeeze(nanmean(nanmean(stress.QRecon.x(virtualGauge.yRange,virtualGauge.xRange,:))));
elseif strcmp('orthotropic',globalOpts.matModel)
    fprintf('WARNING: orthotropic strength identification not implemented.\n')
    
elseif strcmp('orthotropicAngle',globalOpts.matModel)
    fprintf('WARNING: angled orthotropic strength identification not implemented.\n')
    
else
    fprintf('WARNING: specifed material model not recognised.\n')
end

end

