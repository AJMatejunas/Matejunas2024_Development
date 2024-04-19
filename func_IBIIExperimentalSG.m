function SG = func_IBIIExperimentalSG(accel,strain,pos,time,TestDeg,...
    density)
% This function is written to evaluate the stress gage equations for
    %experimentally acquired IBII data. It also produces plots of the
    %results.
    
%Author: Andrew Matejunas

%Data Completed:

%Function Input arguments
    %accel-Structute of measured accelerations with fields
        %x- accelerations in the X direction
        %y- accelerations in the Y direction
    %strain- structure of strains with fields
        %x- normal strains in x direction
        %y- normal strains in y direction
        %s- in-plane shear strains
    %pos- coordinates of specimen with fields
        %x- X coordinates
        %y- Y coordinates
    %Time- Stucture of time data
    %TestDeg- test designation
    %density- density of the material
    
%Output Arguments
    %SG- Stress gauge averaged stresses with fields
        %x- normal stresses in the x direction
        %s- in plane shear stresses

%% Calculate X_vec
X_vec=pos.x;

%% Evaluate stress gage
fprintf('Evaluating Stress Gauge equations \n')
SG=func_Full_SG(accel,X_vec,time,density);

%% Plot Stress Gage Stresses
fprintf('Plotting Stress Gage Stresses \n')

SGplots=func_PlotStressGageStresses(SG,strain,X_vec,time,TestDeg);

% %% Plot Frames for Movie of Stress Wave Propagation
SG_frames=func_GenerateSGMovieFrames(SG,X_vec,time,TestDeg);


end

