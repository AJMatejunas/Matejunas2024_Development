function [specimen,DIC] = func_updateSpecGeomDIC_v4(specimen,material,...
    DIC,disp)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 24/1/2017
% Date Edited: 24/1/2019
% Date Edited: 13/8/2020 - also update the specimen volume and mass based
% on the density
%
%Modified by- Andrew Matejunas (andrew.matejunas@gmail.com)
    %2023-08-02- Converted code to work with DIC

%
% Updates the specimen geometry in the specimen and DIC structs based on
% the selected FOV from the DIC method image processing.

% Check that the nominal specimen dimensions and the dimensions of the
% selected DIC window are the same, if they are not update them

    checkLength = DIC.mPerST*size(disp.x,2);
    if checkLength ~= specimen.length
        specimen.length = checkLength;
        DIC.length = checkLength;
    end
    checkHeight = DIC.mPerST*size(disp.x,1);
    if checkHeight ~= specimen.height
        specimen.height = checkHeight;
        DIC.height = checkHeight;
    end
specimen.volume = specimen.length*specimen.height*specimen.thickness;
specimen.mass = specimen.volume*material.rho;


