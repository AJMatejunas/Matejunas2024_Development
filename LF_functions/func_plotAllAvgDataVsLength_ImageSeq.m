function func_plotAllAvgDataVsLength_ImageSeq(plotParams,imagePath,...
    time,pos,disp,accel,strain,stress)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 9/3/2017
%
% Plot Avg Data vs Length as an Image Sequence

    % Create a plot props structure for nice formatting
    plotProps = func_initPlotPropsStruct();
    
    %----------------------------------------------------------------------
    % Plot Image Sequence of Width Avg Disp X and Strain X
    
    % Create a new directory for this image sequence
    imageSeqSavePath = [imagePath,'ImageSeq_DispXAvg_StrainXAvg\'];
    if exist(imageSeqSavePath,'file') ~= 7
        mkdir(imageSeqSavePath)
    end
    % Assign variables to be plotted
    [deformDisp,rigidDisp] = func_calcDispRemoveRigidBody(disp.x);
    deformDispXAvg = squeeze(mean(deformDisp));
    plotVars{1} = deformDispXAvg*10^3;
    plotVars{2} = strain.xAvg*10^3;
    % Setup the plot parameters and label strings
    plotParams.Rows = 1;
    plotParams.Cols = 2;
    plotParams.axisRange{1} = [0,round(pos.x(end),2),-0.3,0.3];
    plotParams.axisRange{2} = [0,round(pos.x(end),2),-20,20];
    labelStrs.x{1} = 'Length X, L_{s}, (mm)';
    labelStrs.y{1} = 'Displacement X, \delta_{x}, (mm)';
    labelStrs.y{2} = 'Strain X, \epsilon_{x}, (mm/m)';
    % Pass everything to plot the image sequence
    func_plotImageSeqAvgDataVsLength(plotParams,labelStrs,imageSeqSavePath,pos,time,plotVars)
    
    %----------------------------------------------------------------------
    % Plot Image Sequence of Width Avg Disp X and Strain X
    
    % Create a new directory for this image sequence
    imageSeqSavePath = [imagePath,'ImageSeq_AccelXAvg_StressXAvg\'];
    if exist(imageSeqSavePath,'file') ~= 7
        mkdir(imageSeqSavePath)
    end
    % Assign variables to be plotted
    plotVars{1} = accel.xAvg;
    plotVars{2} = stress.xAvg*10^-6;
    % Setup the plot parameters and label strings
    plotParams.Rows = 1;
    plotParams.Cols = 2;
    plotParams.axisRange{1} = [0,round(pos.x(end),2),-7e6,7e6];
    plotParams.axisRange{2} = [0,round(pos.x(end),2),-120,80];
    labelStrs.x{1} = 'Length X, L_{s}, (mm)';
    labelStrs.y{1} = 'Accel X, a_{x}, (m.s^{-2})';
    labelStrs.y{2} = 'Stress X, \sigma_{x}, (MPa)';
    % Pass everything to plot the image sequence
    func_plotImageSeqAvgDataVsLength(plotParams,labelStrs,imageSeqSavePath,pos,time,plotVars)

end