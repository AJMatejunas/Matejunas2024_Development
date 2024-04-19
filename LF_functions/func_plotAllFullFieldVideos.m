function func_plotAllFullFieldVideos(globalOpts,plotParams,imagePath,pos,time,disp,...
    strain,strainRate,accel)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 13/3/2017
%
% Plots 2D frames of data as an image sequence
    
    % Label strings for the image axis
    labelStrs.x = 'X [$mm$]';
    labelStrs.y = 'Y [$mm$]';
    
    %----------------------------------------------------------------------
    % Image Sequence of X Disp, Strain, Accel and Strain Rate

    % Check path for the image sequence exists and if it should be plotted
    imageSeqSavePath = [imagePath,'ImageSeq1_AllFieldsXOnly\'];
    if strcmp(globalOpts.plotImageSeqs,'prompt')
        plotImageSeq = func_checkIfPathExistsDlg(imageSeqSavePath,...
            'X field image sequence folder found, plot again?','Plot Image Sequence?'); 
    elseif strcmp(globalOpts.plotImageSeqs,'yes') || strcmp(globalOpts.plotImageSeqs,'diagnostic')
        plotImageSeq = true;
        if exist(imageSeqSavePath,'file') ~= 7
            mkdir(imageSeqSavePath);
        end
    else
        plotImageSeq = false; 
    end
    
    if plotImageSeq
        % Create title strings and plot variables for the video
        labelStrs.t{1} = 'Displacement $\delta_{x}$, [$mm$]';
        labelStrs.t{2} = 'Accel $a_{x}$, [$m.s^{-2}$]';
        labelStrs.t{3} = 'Strain $\epsilon_{xx}$, [$mm.m^{-1}$]';
        labelStrs.t{4} = 'Strain Rate $\dot{\epsilon_{xx}}$, [$s^{-1}$]';
        % Assign the different variables to be plotted
        [deformDisp.x,~] = func_calcDispRemoveRigidBody(disp.x);
        plotVars{1} = deformDisp.x*10^3;
        plotVars{2} = accel.x;
        plotVars{3} = strain.x*10^3;
        plotVars{4} = strainRate.x;
        % Specify the parameters of the plot
        plotParams.Rows = 2;
        plotParams.Cols = 2;
        plotParams.titleFrameNum = true;
        plotParams.cRange{1} = plotParams.cAxisDisp;
        plotParams.cRange{2} = plotParams.cAxisAccel;
        plotParams.cRange{3} = plotParams.cAxisStrain;
        plotParams.cRange{4} = plotParams.cAxisStrainRate;
        % Plot and save the image sequence to file
        func_plotFullFieldImageSeq(imageSeqSavePath,plotParams,labelStrs,...
            pos,time,plotVars)
        close all
    end
    
    if ~strcmp(globalOpts.plotImageSeqs,'diagnostic')
    %----------------------------------------------------------------------
    % Image Sequence of Disp and Accel - ALL COMPONENTS
    % Create the save path for the heat maps
    plotParams.cutEdgePx = false;
    
    imageSeqSavePath = [imagePath,'ImageSeq2_DispAccelAll\'];
    if strcmp(globalOpts.plotImageSeqs,'prompt')
        plotImageSeq = func_checkIfPathExistsDlg(imageSeqSavePath,...
            'All components of disp/accel image sequence folder found, plot again?','Plot Image Sequence?'); 
    elseif strcmp(globalOpts.plotImageSeqs,'yes')
        plotImageSeq = true;
        if exist(imageSeqSavePath,'file') ~= 7
            mkdir(imageSeqSavePath);
        end
    else
        plotImageSeq = false; 
    end

    if plotImageSeq
        % Calculate the RAW strains for plotting
        [rawStrain.x,~,~] = gradient(disp.x,pos.xStep,pos.yStep,time.step);
        [deformDisp.x,~] = func_calcDispRemoveRigidBody(disp.x);
        [deformDisp.y,~] = func_calcDispRemoveRigidBody(disp.y);
        % Create title strings and plot variables for the video
        labelStrs.t{1} = 'Displacement $\delta_{x}$, [$mm$]';
        labelStrs.t{2} = 'Displacement $\delta_{y}$, [$mm$]';
        labelStrs.t{3} = 'Acceleration $a_{x}$, [$m.s^{-2}$]';
        labelStrs.t{4} = 'Acceleration $a_{y}$, [$m.s^{-2}$]';
        % Assign the different variables to be plotted
        plotVars{1} = deformDisp.x*10^3;
        plotVars{2} = deformDisp.y*10^3;
        plotVars{3} = accel.x;
        plotVars{4} = accel.y;
        % Specify the parameters of the plot
        plotParams.Rows = 2;
        plotParams.Cols = 2;
        plotParams.titleFrameNum = true;
        plotParams.cRange{1} = plotParams.cAxisDisp;
        plotParams.cRange{2} = plotParams.cAxisDisp;
        plotParams.cRange{3} = plotParams.cAxisAccel;
        plotParams.cRange{4} = plotParams.cAxisAccel;
        % Plot and save the image sequence to file
        func_plotFullFieldImageSeq(imageSeqSavePath,plotParams,labelStrs...
            ,pos,time,plotVars)
        close all
    end
    
    %----------------------------------------------------------------------
    % Image Sequence of Strains and Raw Strains - ALL COMPONENTS
    plotParams.cutEdgePx = false;
    
    % Create the save path for the heat maps
    imageSeqSavePath = [imagePath,'ImageSeq3_StrainRaw\'];
    if strcmp(globalOpts.plotImageSeqs,'prompt')
        plotImageSeq = func_checkIfPathExistsDlg(imageSeqSavePath,...
            'All components of strain/raw strain image sequence folder found, plot again?','Plot Image Sequence?');
    elseif strcmp(globalOpts.plotImageSeqs,'yes')
        plotImageSeq = true;
        if exist(imageSeqSavePath,'file') ~= 7
            mkdir(imageSeqSavePath);
        end
    else
        plotImageSeq = false; 
    end

    if plotImageSeq
        % Calculate the RAW strains for plotting
        rawStrain = func_calcStrainFromDisp(disp,pos.xStep,pos.yStep);
        % Create title strings and plot variables for the video
        labelStrs.t{1} = 'Strain $\epsilon_{xx}$, [$mm.m^{-1}$]';
        labelStrs.t{2} = 'Strain $\epsilon_{yy}$, [$mm.m^{-1}$]';
        labelStrs.t{3} = 'Strain $\epsilon_{xy}$, [$mm.m^{-1}$]';
        labelStrs.t{4} = 'Raw Strain $\epsilon_{xx}$, [$mm.m^{-1}$]';
        labelStrs.t{5} = 'Raw Strain $\epsilon_{yy}$, [$mm.m^{-1}$]';
        labelStrs.t{6} = 'Raw Strain $\epsilon_{xy}$, [$mm.m^{-1}$]';
        % Assign the different variables to be plotted
        plotVars{1} = strain.x*10^3;
        plotVars{2} = strain.y*10^3;
        plotVars{3} = strain.s*10^3;
        plotVars{4} = rawStrain.x*10^3;
        plotVars{5} = rawStrain.y*10^3;
        plotVars{6} = rawStrain.s*10^3;
        % Specify the parameters of the plot
        plotParams.Rows = 2;
        plotParams.Cols = 3;
        plotParams.titleFrameNum = true;
        plotParams.cRange{1} = plotParams.cAxisStrain;
        plotParams.cRange{2} = plotParams.cAxisStrain;
        plotParams.cRange{3} = plotParams.cAxisStrain;
        plotParams.cRange{4} = plotParams.cAxisRawStrain;
        plotParams.cRange{5} = plotParams.cAxisRawStrain;
        plotParams.cRange{6} = plotParams.cAxisRawStrain;
        % Plot and save the image sequence to file
        func_plotFullFieldImageSeq(imageSeqSavePath,plotParams,labelStrs,...
            pos,time,plotVars)
        close all
    end 
    
    if strcmp(globalOpts.fieldComponents,'all')
        %----------------------------------------------------------------------
        % Image Sequence of Strains and Strain Rate - ALL COMPONENTS
        plotParams.cutEdgePx = true;

        % Create the save path for the heat maps
        imageSeqSavePath = [imagePath,'ImageSeq4_StrainStrainRate\'];
        if strcmp(globalOpts.plotImageSeqs,'prompt')
            plotImageSeq = func_checkIfPathExistsDlg(imageSeqSavePath,...
                'All components of strain/strain rate image sequence folder found, plot again?','Plot Image Sequence?');   
        elseif strcmp(globalOpts.plotImageSeqs,'yes')
            plotImageSeq = true;
            if exist(imageSeqSavePath,'file') ~= 7
                mkdir(imageSeqSavePath);
            end
        else
            plotImageSeq = false; 
        end
        
        if plotImageSeq   
            % Create title strings and plot variables for the video
            labelStrs.t{1} = 'Strain $\epsilon_{xx}$, [$mm.m^{-1}$]';
            labelStrs.t{2} = 'Strain $\epsilon_{yy}$, [$mm.m^{-1}$]';
            labelStrs.t{3} = 'Strain $\epsilon_{xy}$, [$mm.m^{-1}$]';
            labelStrs.t{4} = 'Strain Rate $\dot{\epsilon_{xx}}$, [$s^{-1}$]';
            labelStrs.t{5} = 'Strain Rate $\dot{\epsilon_{yy}}$, [$s^{-1}$]';
            labelStrs.t{6} = 'Strain Rate $\dot{\epsilon_{xy}}$, [$s^{-1}$]';
            % Assign the different variables to be plotted
            plotVars{1} = strain.x*10^3;
            plotVars{2} = strain.y*10^3;
            plotVars{3} = strain.s*10^3;
            plotVars{4} = strainRate.x;
            plotVars{5} = strainRate.y;
            plotVars{6} = strainRate.s;
            % Specify the parameters of the plot
            plotParams.Rows = 2;
            plotParams.Cols = 3;
            plotParams.titleFrameNum = true;
            plotParams.cRange{1} = plotParams.cAxisStrain;
            plotParams.cRange{2} = plotParams.cAxisStrain;
            plotParams.cRange{3} = plotParams.cAxisStrain;
            plotParams.cRange{4} = plotParams.cAxisStrainRate;
            plotParams.cRange{5} = plotParams.cAxisStrainRate;
            plotParams.cRange{6} = plotParams.cAxisStrainRate;
            % Plot and save the image sequence to file
            func_plotFullFieldImageSeq(imageSeqSavePath,plotParams,labelStrs,...
                pos,time,plotVars)
            close all
        end 
    end
    end

end
