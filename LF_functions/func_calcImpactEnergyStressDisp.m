function EI = func_calcImpactEnergyStressDisp(EBalOpts,pos,specimen,disp,stress,t)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 26/8/2017
% Calculates the impact energy components on the 

    % Check if we have a 3D stress matrix if not then t is always 1
    if length(size(stress)) < 3
        st = 1;
    else
        st = t;
    end
    
    % Remove the options from the struct to shorten code length
    sv = EBalOpts.sv;
    xSubset = EBalOpts.xSubset;
    xBound = EBalOpts.xBound;
    ySubset = EBalOpts.ySubset;
    yBound = EBalOpts.yBound;
    
    if t == 1
        % Calculate the impact energy from the stress
        % LHS 
        EI(1) = sv(1)*sum(sum((stress.x(ySubset,xBound(1),st)*pos.yStep*specimen.thickness)...
            .*disp.x(ySubset,xBound(1),t)));        % xx, LS
        EI(2) = sv(2)*sum(sum((stress.s(ySubset,xBound(1),st)*pos.yStep*specimen.thickness)...
            .*disp.y(ySubset,xBound(1),t)));        % xy, LS
        % RHS 
        EI(3) = sv(3)*sum(sum((stress.x(ySubset,xBound(2),st)*pos.yStep*specimen.thickness)...
            .*disp.x(ySubset,xBound(2),t)));     % xx, RS
        EI(4) = sv(4)*sum(sum((stress.s(ySubset,xBound(2),st)*pos.yStep*specimen.thickness)...
            .*disp.y(ySubset,xBound(2),t)));    % xy, RS 
        % Top
        EI(5) = sv(5)*sum(sum((stress.y(yBound(1),xSubset,st)*pos.xStep*specimen.thickness)...
            .*disp.y(yBound(1),xSubset,t)));         % yy, TS
        EI(6) = sv(6)*sum(sum((stress.s(yBound(1),xSubset,st)*pos.xStep*specimen.thickness)...
            .*disp.x(yBound(1),xSubset,t)));        % yx, TS
        % Bottom - only count if using both data halves
        EI(7) = sv(7)*sum(sum((stress.y(yBound(2),xSubset,st)*pos.xStep*specimen.thickness)...
            .*disp.y(yBound(2),xSubset,t)));     % yy, BS
        EI(8) = sv(8)*sum(sum((stress.s(yBound(2),xSubset,st)*pos.xStep*specimen.thickness)...
            .*disp.x(yBound(2),xSubset,t)));    % yx, BS
    else
        % LHS
        EI(1) = sv(1)*sum(sum(stress.x(ySubset,xBound(1),st)*pos.yStep*specimen.thickness...
            .*(disp.x(ySubset,xBound(1),t)-disp.x(ySubset,xBound(1),t-1))));        % xx, LS
        EI(2) = sv(2)*sum(sum(stress.s(ySubset,xBound(1),st)*pos.yStep*specimen.thickness...
            .*(disp.y(ySubset,xBound(1),t)-disp.y(ySubset,xBound(1),t-1))));        % xy, LS
        % RHS
        EI(3) = sv(3)*sum(sum(stress.x(ySubset,xBound(2),st)*pos.yStep*specimen.thickness...
            .*(disp.x(ySubset,xBound(2),t)-disp.x(ySubset,xBound(2),t-1))));        % xx, RS
        EI(4) = sv(4)*sum(sum(stress.s(ySubset,xBound(2),st)*pos.yStep*specimen.thickness...
            .*(disp.y(ySubset,xBound(2),t)-disp.y(ySubset,xBound(2),t-1))));        % xy, RS
        % Top
        EI(5) = sv(5)*sum(sum(stress.y(yBound(1),xSubset,st)*pos.xStep*specimen.thickness...
            .*(disp.y(yBound(1),xSubset,t)-disp.y(yBound(1),xSubset,t-1))));        % yy, TS
        EI(6) = sv(6)*sum(sum(stress.s(yBound(1),xSubset,st)*pos.xStep*specimen.thickness...
            .*(disp.x(yBound(1),xSubset,t)-disp.x(yBound(1),xSubset,t-1))));        % yx, TS
        % Bottom - only count if using both data halves
        EI(7) = sv(7)*sum(sum(stress.y(yBound(2),xSubset,st)*pos.xStep*specimen.thickness...
            .*(disp.y(yBound(2),xSubset,t)-disp.y(yBound(2),xSubset,t-1))));        % yy, BS
        EI(8) = sv(8)*sum(sum(stress.s(yBound(2),xSubset,st)*pos.xStep*specimen.thickness...
            .*(disp.x(yBound(2),xSubset,t)-disp.x(yBound(2),xSubset,t-1))));        % yx, BS              
    end
end

