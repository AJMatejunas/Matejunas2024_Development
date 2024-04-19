function stress = func_calcStressFromStrain2D_SingleFrame(strain,Q,t)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 26/8/2017
% Calculates the 2D stress components using the stiffness matrix
    
    [sy,sx,~] = size(strain.x);
    for x = 1:sx
        for y = 1:sy
            stress.x(y,x) = Q(1,1)*strain.x(y,x,t)+...
                Q(1,2)*strain.y(y,x,t)+Q(1,3)*strain.s(y,x,t);

            stress.y(y,x) = Q(2,1)*strain.x(y,x,t)+...
                Q(2,2)*strain.y(y,x,t)+Q(2,3)*strain.s(y,x,t);

            stress.s(y,x) = Q(3,1)*strain.x(y,x,t)+...
                Q(3,2)*strain.y(y,x,t)+Q(3,3)*strain.s(y,x,t);
        end
    end
end

