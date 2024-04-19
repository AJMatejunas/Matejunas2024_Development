function func_plotAllFullFieldVideos_EBal(plotParams,imagePath,pos,time,disp,...
    vel,accel,strain)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 13/3/2017
%
% Plots 2D frames of data as an image sequence
    
    plotParams.cAxisType = 'MaxQuantile';
    % Label strings for the image axis
    labelStrs.x = 'X (mm)';
    labelStrs.y = 'Y (mm)';

    %----------------------------------------------------------------------
    % Image Sequence of X,Y Disp and Vel
    
    % Create the save path for the heat maps
    imageSeqSavePath = [imagePath,'ImageSeq_DispVelXY\'];

    % Check if the image sequence path exits, if it does ask the user if it
    % needs to be plotted again
    if exist(imageSeqSavePath,'file') == 7
        choice = questdlg('Full-field disp,vel XY sequence folder found, plot again?', ...
        'Plot Image Sequence?', 'Yes','No','No');
        switch choice
            case 'Yes'
                plotImageSeq = true;
            case 'No'
                plotImageSeq = false;
        end
    else
        % If the directory does not exist, create it and plot the image seq
        mkdir(imageSeqSavePath);
        plotImageSeq = true;
    end
    
    if plotImageSeq
        % Create title strings and plot variables for the video
        labelStrs.t{1} = 'Displacement X, \delta_{x}, (mm)';
        labelStrs.t{3} = 'Displacement Y, \delta_{y}, (mm)';
        labelStrs.t{2} = 'Velocity X, v_{x}, (m.s^{-1})';
        labelStrs.t{4} = 'Velocity Y, v_{y}, (m.s^{-1})';
        % Assign the different variables to be plotted
        [deformDisp,rigidDisp] = func_calcDispRemoveRigidBody(disp.x);
        plotVars{1} = deformDisp*10^3;
        plotVars{3} = disp.y*10^3;
        plotVars{2} = vel.x;
        plotVars{4} = vel.y;
        % Specify the parameters of the plot
        plotParams.Rows = 2;
        plotParams.Cols = 2;
        plotParams.titleFrameNum = true;
        plotParams.cRange{1} = plotParams.cAxisDisp;
        plotParams.cRange{3} = plotParams.cAxisDisp;
        plotParams.cRange{2} = plotParams.cAxisVel;
        plotParams.cRange{4} = plotParams.cAxisVel;
        % Plot and save the image sequence to file
        func_plotFullFieldImageSeq(imageSeqSavePath,plotParams,labelStrs,pos,time,plotVars) 
    end
    
    %----------------------------------------------------------------------
    % Image Sequence of X,Y,S Strains
    
    % Create the save path for the heat maps
    imageSeqSavePath = [imagePath,'ImageSeq_StrainXY\'];

    % Check if the image sequence path exits, if it does ask the user if it
    % needs to be plotted again
    if exist(imageSeqSavePath,'file') == 7
        choice = questdlg('Full-field strain X,Y,XY image sequence folder found, plot again?', ...
        'Plot Image Sequence?', 'Yes','No','No');
        switch choice
            case 'Yes'
                plotImageSeq = true;
            case 'No'
                plotImageSeq = false;
        end
    else
        % If the directory does not exist, create it and plot the image seq
        mkdir(imageSeqSavePath);
        plotImageSeq = true;
    end
    
    if plotImageSeq
        % Create title strings and plot variables for the video
        labelStrs.t{1} = 'Strain X, \epsilon_{x}, (mm.m^{-1})';
        labelStrs.t{2} = 'Strain Y, \epsilon_{y}, (mm.m^{-1})';
        labelStrs.t{3} = 'Strain XY, \epsilon_{xy}, (mm.m^{-1})';
        % Assign the different variables to be plotted
        plotVars{1} = strain.x*10^3;
        plotVars{2} = strain.y*10^3;
        plotVars{3} = strain.s*10^3;
        % Specify the parameters of the plot
        plotParams.Rows = 3;
        plotParams.Cols = 1;
        plotParams.titleFrameNum = true;
        % Plot and save the image sequence to file
        func_plotFullFieldImageSeq(imageSeqSavePath,plotParams,labelStrs,pos,time,plotVars) 
    end
end