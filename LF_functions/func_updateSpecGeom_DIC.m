function [specimen,DIC] = func_updateSpecGeom_DIC(specimen,DIC,disp)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 24/1/2017
% Date Edited: 24/1/2019
%Modified by- Andrew Matejunas (andrew.matejunas@gmail.com)
    %2023-08-02- Converted code to work with DIC

%
% Updates the specimen geometry in the specimen and DIC structs based on
% the selected DIC from the MatchID DIC processing.

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




