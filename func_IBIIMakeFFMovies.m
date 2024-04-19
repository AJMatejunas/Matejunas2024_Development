function func_IBIIMakeFFMovies(pos,time,disp,...
    SmoothStrain,SmoothAccel,smoothOpts,extrapOpts,diffOpts,ExpDesig,SavePath)
%This script is written to make full field movies of displacement, strain,
    %and accelerations for an IBII test. 

%Author: Andrew Matejunas

%Version History/Change log:
    %2023-06-24: Original version

%Function input arguments
    %pos- structure containing coordinates
    %time- structure containing time information
    %disp- structure containing displacement components with and without
        %smoothing
    %SmoothStrain- structure containing smoothed and corrected strain
        %fields
    %SmoothAccel- Structure containing smoothed and corrected acceleration
        %fields
    %SmoothOpts- Structure containing smoothing information
    %extrapOpts- Structure containing kinematic field edge extrapolation
        %options
    %diffOpts- structure containing temporal differentiation options for
        %calculation of accelerations
    %ExpDesig- Test designation
    %SavePath- Path to save movie frames

%Function outputs
    %Frames for a Full-field movie of Displacements
    %Frames for a movie of Full-field strains
    %Frames for a movie of full field accelerations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create Directory for the movies
MovieDir=strcat(SavePath,'/Vids');
mkdir(MovieDir)

%% Create displacement video
func_IBIIplotFFdisp(pos,time,disp,smoothOpts,ExpDesig,MovieDir);

%% Create Strain video
func_IBIIplotFFstrainMovie(pos,time,disp,SmoothStrain,smoothOpts,...
    extrapOpts,ExpDesig,MovieDir);

%% Create Acceleration Video
func_IBIIplotFFaccelMovie(pos,time,disp,SmoothAccel,smoothOpts,extrapOpts,...
    diffOpts,ExpDesig,MovieDir)

end