function [PPTPlot,PubPlot,FSPlot] =...
    func_PlotIdentErrorContours(SpaKernVec,TempKernVec,IdentErrors,...
    ColorScheme,PPTsize,PubSize,...
    SaveDir,ParentDesig)

%This function is written to produce identification error heat maps for
    %the dynamic part of a viscoelastic material characterized by the
    %Maxwell formulation of the standard solid model. This version produces
    %filled contour plots

%Author: Andrew Matejunas
%Date Completed: 2022-12-02
%Version history/Change log:

%Function Input Arguments
    %SpaKernVec- Vector of spatial smoothing kernal sizes (pixels)
    %TempKernVec- Vector of temporal smoothing kernal (sizes) frames
    %IdentErrors- structure containing constitutive parameter
        %identification errors (%) with fields 
            %K- Bulk modulus errors [# Temporal Kernals, # Spatial Kernals]
            %G- Shear modulus identification errors 
                % [# Temporal Kernals, # Spatial Kernals]
            %tau- Time constant identification errors 
                % [# Temporal Kernals, # Spatial Kernals]
            %RMS- Root mean square errors 
                % [# Temporal Kernals, # Spatial Kernals]
   %ColorScheme- string defining the color map for the contour plot
   %PPT size- size of the power point version of the figure [

%Function Output Arguments:


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outputArg1 = inputArg1;
outputArg2 = inputArg2;
end