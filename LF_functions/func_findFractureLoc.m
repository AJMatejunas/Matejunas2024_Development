function fracture = func_findFractureLoc(imageSeqPath,pos,time,disp,accel,strain,strainRate)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 27/3/2017
% Date Edited: 11/4/2019
%
% Prompts user to select frames in which crack growth occurs, asks user to
% select the location of the crack on the displacement or strain maps.
% Returns all information in the 'fracture' struct
    
    % Prompt the user with the location of the video of disp / strain
    msg = {'Image sequences of kinematic fields are located:',imageSeqPath,...
        '','Find fracture frame using raw strain maps.'};
    h_mdlg = msgbox(msg,'Image Location');
    uiwait(h_mdlg);
    
    % Ask the user for the frame to select the crack location 
    while true
        frameData = inputdlg({'Frame to select fracture location:'}, ...
                 'Fracture Frame Selection', 1, {num2str(time.numFrames)} );
        fracture.locFrame = str2double(frameData{1});
        % Make sure the user is giving a sensible input
        if (fracture.locFrame<=time.numFrames) && (fracture.locFrame>=1)
            break;
        else
            msg = 'Selected frame does not exist!';
            h_mdlg = msgbox(msg,'Warning!');
            uiwait(h_mdlg);
        end
    end
    
    % Create the plot properties structure
    plotProps = func_initPlotPropsStruct();
    
    % Calculate the RAW strains for plotting
    xStep = abs(pos.x(2) - pos.x(1));
    yStep = abs(pos.y(2) - pos.y(1));
    [rawStrain.x,~,~] = gradient(disp.x,xStep,yStep,time.step);
    
    % Number of plots and rows/cols for subplot
    numPlots = 4;
    rows = 2;
    cols = 2;
    % Get the deformation portion of the displacement
    [deformDisp,~] = func_calcDispRemoveRigidBody(disp.x);
    % Setup the variables for plotting the heat maps
    mapVar{1} = deformDisp*10^3;
    mapVar{2} = accel.x;
    mapVar{3} = rawStrain.x*10^3;
    mapVar{4} = strain.x*10^3;
    %mapVar{4} = strainRate.x;
    % Title strings for each figure
    titleStr{1} = {'Select crack origin location:','Displacement, $\delta_{x}$, ($mm$)'};
    titleStr{2} = {'Select crack origin location:','Accel, $a_{x}$, ($m.s^{-2}$)'};
    titleStr{3} = {'Select crack origin location:','Raw Strain, $\epsilon_{xx}$, ($mm.m^{-1}$)'};
    titleStr{4} = {'Select crack origin location:','Strain, $\epsilon_{xx}$, ($mm.m^{-1}$)'};
    %titleStr{4} = {'Select crack origin location:','Strain Rate X, d\epsilon_{xx}/dt, (m^{-1})'};
    % Find the caxis range for the displacement and strain
    cRange{1} = func_calcColourBarRange('MaxQuantile',deformDisp,[0.1,0.9],1)*10^3;
    cRange{2} = func_calcColourBarRange('MaxQuantile',accel.x,[0.05,0.95],1.2);
    cRange{3} = func_calcColourBarRange('MaxQuantile',rawStrain.x,[0.05,0.95],1.2)*10^3;
    cRange{4} = func_calcColourBarRange('MaxQuantile',strain.x,[0.05,0.95],1.2)*10^3;
    %cRange{4} = func_calcColourBarRange('MaxQuantile',strainRate.x,[0.05,0.95],1.2);
     
    % Create and size the figure
    hf = figure('Visible','On');
    set(hf,'Position', [plotProps.screenSize(1),plotProps.screenSize(2),...
        plotProps.screenSize(4)*1.6,plotProps.screenSize(4)])
    set(hf,'PaperPositionMode','auto')
    % Loop until crack location is accepted
    while true
        % Plot the heat maps for selecting the crack 
        for p = 1:numPlots
            hsp{p} = subplot(rows,cols,p);
            imagesc(squeeze(mapVar{p}(:,:,fracture.locFrame)))
            colorbar
            colormap(jet)
            caxis(cRange{p}); 
            set(gca,'fontsize', plotProps.fs,'fontname',plotProps.ft)
            title(titleStr{p},'fontsize',plotProps.hfs,'fontname',plotProps.ft,'Interpreter',plotProps.interpText)
            axis image
        end
        
        % Tell the user to select the crack location
        msg = {'Select crack origin location.'};
        h_mdlg = msgbox(msg,'Image Location');
        uiwait(h_mdlg);
        
        % Use ginput to get the crack location
        [xLoc,yLoc] = ginput(1);
        
        % Use the plot command to show the selected location on the figure
        for p = 1:numPlots
            subplot(hsp{p})
            hold on
            plot(xLoc,yLoc,'+k','markersize',15)
            title(sprintf('X = %ipx, Y = %ipx',round(xLoc),round(yLoc)))
            hold off
        end
        
        % Prompt to accept or repeat
        choice = questdlg('Is the selected crack location ok?',...
            'Process Grid Images','Yes','No','Yes');
        switch choice
            case 'Yes'
                fracture.locX = round(xLoc);
                fracture.locY = round(yLoc);
                break
            case 'No'
                clf(hf);
                continue
        end
    end
    close(hf);  
end