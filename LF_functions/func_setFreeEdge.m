function disp = func_setFreeEdge(freeEdge,disp)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 30/8/2017
%
% Flips displacement fields based on specimen orientation

    if strcmp(freeEdge,'Right')
        % If the free edge is on the right hand side flip the displacement
        % matrix left to right
        fprintf('Free edge is on the right hand side, flipping displacement matrices\n')
        for f = 1:size(disp.x,3)
            disp.x(:,:,f) = -fliplr(disp.x(:,:,f));
            disp.y(:,:,f) = fliplr(disp.y(:,:,f));
        end
    elseif strcmp(freeEdge,'Left')
        fprintf('Free edge is on the left hand side, no correction required.\n')
    end

end