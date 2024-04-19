function [stressX,stressY,stressS] = func_calcStressFromStrain2D(Q,strainX,strainY,strainS)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 8/8/2017
% Date Modified: 1/9/2018
% Calculates the 2D stress components using the stiffness matrix
    
    [sy,sx,st] = size(strainX);
    stressX = zeros(sy,sx,st);
    stressY = zeros(sy,sx,st);
    stressS = zeros(sy,sx,st);
    for t = 1:st
        for x = 1:sx
            for y = 1:sy
                stressX(y,x,t) = Q(1,1)*strainX(y,x,t)+...
                    Q(1,2)*strainY(y,x,t)+Q(1,3)*strainS(y,x,t);

                stressY(y,x,t) = Q(2,1)*strainX(y,x,t)+...
                    Q(2,2)*strainY(y,x,t)+Q(2,3)*strainS(y,x,t);

                stressS(y,x,t) = Q(3,1)*strainX(y,x,t)+...
                    Q(3,2)*strainY(y,x,t)+Q(3,3)*strainS(y,x,t);
            end
        end
    end
end

