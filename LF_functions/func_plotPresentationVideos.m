function func_plotPresentationVideos(globalOpts,plotParams,savePath,pos,time,...
    strain,strainRate,accel,stress,identStiff,ssCurve)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 29/3/2017
%
% Plots image sequences used to create videos for power points presentations 

% Create the save paths for the image sequences
imageSeqSavePath{1} = [savePath,'Presentation_ImageSeq1_MapsSSCurve\'];
%imageSeqSavePath{2} = [savePath,'Presentation_ImageSeq2_MapsSSCurveFit\'];
for p = 1:length(imageSeqSavePath)
    if exist(imageSeqSavePath{p},'file') ~= 7
        mkdir(imageSeqSavePath{p});
    end
end
plotInd = 0;

%--------------------------------------------------------------------------
% Depending on which material model we use plot the SS curve      
if strcmp('orthotropicReduced',globalOpts.matModel)
    % Video 1: Strain Map and Stress Strain Curve
    plotInd = 1; 
    plotParams.Rows = 2;
    plotParams.Cols = 2;
    % Options for the stress/strain curve being plotted
    ssCurve.stress = stress.xAvg(ssCurve.locX,:);
    ssCurve.strain = strain.xAvg(ssCurve.locX,:);
    ssCurve.linFitCoeffs = identStiff.QxxLinFitCoeffs{ssCurve.locX};
    ssCurve.plotTrendLine = false; 
    % Axis label strings for each plot
    labelStrs.x{1} = '$\overline{\epsilon_{xx}}^{y}$ [$mm.m^{-1}$]';
    labelStrs.x{2} = 'X [$mm$]';
    labelStrs.x{3} = 'X [$mm$]';
    labelStrs.x{4} = 'X [$mm$]';
    labelStrs.y{1} = '$\overline{\sigma_{xx}}^{y}$ [$MPa$]';
    labelStrs.y{2} = 'Y [$mm$]';
    labelStrs.y{3} = 'Y [$mm$]';
    labelStrs.y{4} = 'Y [$mm$]';
    % Create title strings and plot variables for the video
    labelStrs.t{1} = ['(a) x = ',sprintf('%0.2f',pos.x(ssCurve.locX)*10^3),' $mm$'];
    labelStrs.t{2} = '(b) $a_{x}$ [$m.s^{-2}$]';
    labelStrs.t{3} = '(c) $\epsilon_{xx}$ [$mm.m^{-1}$]';
    labelStrs.t{4} = '(d) $\dot{\epsilon_{xx}}$ [$s^{-1}$]';
    % Assign the different map variables to be plotted
    mapVars{2} = accel.x;
    mapVars{3} = strain.x*10^3;
    mapVars{4} = strainRate.x;
    % Specify colour bar ranges for the plot
    plotParams.cRange{2} = plotParams.cAxisAccel;
    plotParams.cRange{3} = plotParams.cAxisStrain;
    plotParams.cRange{4} = plotParams.cAxisStrainRate;
    
    fprintf('\tPlotting stress-strain curve and variable maps for full time range.\n')
    func_plotMapsAndSSCurveGeneric(imageSeqSavePath{plotInd},plotParams,labelStrs,...
        pos,time,ssCurve,mapVars);
    
    % Video 2: Strain Map and Stress Strain Curve
    %{
    fprintf('\tPlotting stress-strain curve and variable maps for compression loading only.\n')
    plotInd = 2; 
    ssCurve.plotTrendLine = true;
    [~,minFrame] = min(ssCurve.stress); 
    plotParams.tRange = 1:minFrame;
    func_plotMapsAndSSCurveGeneric(imageSeqSavePath{plotInd},plotParams,labelStrs,...
        pos,time,ssCurve,mapVars);
    %}
    
elseif strcmp('orthotropic',globalOpts.matModel) || strcmp('isotropic',globalOpts.matModel)
    % Video 1: Strain Map and Stress Strain Curve
    plotInd = 1; 
    plotParams.Rows = 2;
    plotParams.Cols = 2;
    % Options for the stress/strain curve being plotted
    ssCurve.stress = stress.xAvg(ssCurve.locX,:);
    ssCurve.strain = strain.xnyAvg(ssCurve.locX,:);
    ssCurve.linFitCoeffs = identStiff.QxxLinFitCoeffs{ssCurve.locX};
    ssCurve.plotTrendLine = false; 
    % Axis label strings for each plot
    labelStrs.x{1} = '$\overline{\epsilon_{xx}+\nu_{xy} \epsilon_{yy}}^{y}$ [$mm.m^{-1}$]';
    labelStrs.x{2} = 'X [$mm$]';
    labelStrs.x{3} = 'X [$mm$]';
    labelStrs.x{4} = 'X [$mm$]';
    labelStrs.y{1} = '$\overline{\sigma_{xx}}^{y}$ [$MPa$]';
    labelStrs.y{2} = 'Y [$mm$]';
    labelStrs.y{3} = 'Y [$mm$]';
    labelStrs.y{4} = 'Y [$mm$]';
    % Create title strings and plot variables for the video
    labelStrs.t{1} = ['(a) x = ',sprintf('%0.2f',pos.x(ssCurve.locX)*10^3),' $mm$'];
    labelStrs.t{2} = '(b) $a_{x}$ [$m.s^{-2}$]';
    labelStrs.t{3} = '(c) $\epsilon_{xx}$ [$mm.m^{-1}$]';
    labelStrs.t{4} = '(d) $\dot{\epsilon_{xx}}$ [$s^{-1}$]';
    % Assign the different map variables to be plotted
    mapVars{2} = accel.x;
    mapVars{3} = strain.x*10^3;
    mapVars{4} = strainRate.x;
    % Specify colour bar ranges for the plot
    plotParams.cRange{2} = plotParams.cAxisAccel;
    plotParams.cRange{3} = plotParams.cAxisStrain;
    plotParams.cRange{4} = plotParams.cAxisStrainRate;
    
    fprintf('\tPlotting stress-strain curve and variable maps for full time range.\n')
    func_plotMapsAndSSCurveGeneric(imageSeqSavePath{plotInd},plotParams,labelStrs,...
        pos,time,ssCurve,mapVars);
    
    % Video 2: Strain Map and Stress Strain Curve
    %{
    fprintf('\tPlotting stress-strain curve and variable maps for compression loading only.\n')
    plotInd = 2; 
    ssCurve.plotTrendLine = true;
    [~,minFrame] = min(ssCurve.stress); 
    plotParams.tRange = 1:minFrame;
    func_plotMapsAndSSCurveGeneric(imageSeqSavePath{plotInd},plotParams,labelStrs,...
        pos,time,ssCurve,mapVars);
    %}
    
elseif strcmp('orthotropicAngle',globalOpts.matModel)
    fprintf('WARNING: angled orthotropic strength identification not implemented.\n')
else
    fprintf('WARNING: specifed material model not recognised.\n')
end 

end

