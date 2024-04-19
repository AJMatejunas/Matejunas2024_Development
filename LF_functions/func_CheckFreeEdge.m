function [freeEdge,specimen,disp] = func_CheckFreeEdge(...
    specimen,disp,printToCons)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 8/12/2017
% Date Edited: 3/5/2019

%Edited by Andrew Matejunas 08/03/2023
    %Reduced complexity of the code to work with DIC data obtained from
    %Match ID. Removed ability to select a free edge Now the free edge is
    %hardcoded and flipped if needed.
%
% Determines which side of the reference image contains the free edge of
% the specimen for stress gauge processing
        
    % Allow printing to the console by default
    if nargin < 3
        printToCons = true;
    end

freeEdge=specimen.freeEdge;

    if strcmp(freeEdge,'Right')
        % If the free edge is on the right hand side flip the displacement
        % matrix left to right
        if printToCons
            fprintf('Free edge is on the right hand side, flipping displacement matrices\n')
        end
        for f = 1:size(disp.x,3)
            disp.x(:,:,f) = -fliplr(disp.x(:,:,f));
            disp.y(:,:,f) = fliplr(disp.y(:,:,f));
        end
    elseif strcmp(freeEdge,'Left')
        if printToCons
            fprintf('Free edge is on the left hand side, no correction required.\n')
        end
    end
end