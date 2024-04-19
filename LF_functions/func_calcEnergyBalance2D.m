function energy = func_calcEnergyBalance2D(EBalOpts,pos,time,specimen,element...
    ,disp,vel,accel,strain,stress)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 8/8/2017
% Calculates the kinetic, strain and impact energy components using input
% variables.
    
    % Remove the options from the struct to shorten code length
    sv = EBalOpts.sv;
    xSubset = EBalOpts.xSubset;
    xBound = EBalOpts.xBound;
    ySubset = EBalOpts.ySubset;
    yBound = EBalOpts.yBound;
    impEdgeX = EBalOpts.impEdgeX;
    impEdgeY = EBalOpts.impEdgeY;
    
    for t = 1:time.numFrames
        %------------------------------------------------------------------
        % EXTERNAL WORK - Calculate from const law 
        if t == 1
            % Calculate the impact energy from the stress
            % LHS 
            EI(1) = sv(1)*sum(sum((stress.x(ySubset,xBound(1),t)*pos.yStep*specimen.thickness)...
                .*disp.x(ySubset,xBound(1),t)));        % xx, LS
            EI(2) = sv(2)*sum(sum((stress.s(ySubset,xBound(1),t)*pos.yStep*specimen.thickness)...
                .*disp.y(ySubset,xBound(1),t)));        % xy, LS
            % RHS 
            EI(3) = sv(3)*sum(sum((stress.x(ySubset,xBound(2),t)*pos.yStep*specimen.thickness)...
                .*disp.x(ySubset,xBound(2),t)));     % xx, RS
            EI(4) = sv(4)*sum(sum((stress.s(ySubset,xBound(2),t)*pos.yStep*specimen.thickness)...
                .*disp.y(ySubset,xBound(2),t)));    % xy, RS 
            % Top
            EI(5) = sv(5)*sum(sum((stress.y(yBound(1),xSubset,t)*pos.xStep*specimen.thickness)...
                .*disp.y(yBound(1),xSubset,t)));         % yy, TS
            EI(6) = sv(6)*sum(sum((stress.s(yBound(1),xSubset,t)*pos.xStep*specimen.thickness)...
                .*disp.x(yBound(1),xSubset,t)));        % yx, TS
            % Bottom - only count if using both data halves
            EI(7) = sv(7)*sum(sum((stress.y(yBound(2),xSubset,t)*pos.xStep*specimen.thickness)...
                .*disp.y(yBound(2),xSubset,t)));     % yy, BS
            EI(8) = sv(8)*sum(sum((stress.s(yBound(2),xSubset,t)*pos.xStep*specimen.thickness)...
                .*disp.x(yBound(2),xSubset,t)));    % yx, BS

            energy.impactFromStrain(t) = sum(EI);
            clear EI      
        else
            % LHS
            EI(1) = sv(1)*sum(sum(stress.x(ySubset,xBound(1),t)*pos.yStep*specimen.thickness...
                .*(disp.x(ySubset,xBound(1),t)-disp.x(ySubset,xBound(1),t-1))));        % xx, LS
            EI(2) = sv(2)*sum(sum(stress.s(ySubset,xBound(1),t)*pos.yStep*specimen.thickness...
                .*(disp.y(ySubset,xBound(1),t)-disp.y(ySubset,xBound(1),t-1))));        % xy, LS
            % RHS
            EI(3) = sv(3)*sum(sum(stress.x(ySubset,xBound(2),t)*pos.yStep*specimen.thickness...
                .*(disp.x(ySubset,xBound(2),t)-disp.x(ySubset,xBound(2),t-1))));        % xx, RS
            EI(4) = sv(4)*sum(sum(stress.s(ySubset,xBound(2),t)*pos.yStep*specimen.thickness...
                .*(disp.y(ySubset,xBound(2),t)-disp.y(ySubset,xBound(2),t-1))));        % xy, RS
            % Top
            EI(5) = sv(5)*sum(sum(stress.y(yBound(1),xSubset,t)*pos.xStep*specimen.thickness...
                .*(disp.y(yBound(1),xSubset,t)-disp.y(yBound(1),xSubset,t-1))));        % yy, TS
            EI(6) = sv(6)*sum(sum(stress.s(yBound(1),xSubset,t)*pos.xStep*specimen.thickness...
                .*(disp.x(yBound(1),xSubset,t)-disp.x(yBound(1),xSubset,t-1))));        % yx, TS
            % Bottom - only count if using both data halves
            EI(7) = sv(7)*sum(sum(stress.y(yBound(2),xSubset,t)*pos.xStep*specimen.thickness...
                .*(disp.y(yBound(2),xSubset,t)-disp.y(yBound(2),xSubset,t-1))));        % yy, BS
            EI(8) = sv(8)*sum(sum(stress.s(yBound(2),xSubset,t)*pos.xStep*specimen.thickness...
                .*(disp.x(yBound(2),xSubset,t)-disp.x(yBound(2),xSubset,t-1))));        % yx, BS              

            if EBalOpts.impactEnerCutOff && (t>EBalOpts.impactEnerCutOffFrame)
                energy.impactFromStrain(t) = energy.impactFromStrain(t-1);
            else
                energy.impactFromStrain(t) = energy.impactFromStrain(t-1)+sum(EI);
            end
            clear EI
        end
        
        %------------------------------------------------------------------
        % EXTERNAL WORK - Calculate from Acceleration
        if EBalOpts.calcExtWorkFromAccel
            energy.forceX(t) = specimen.mass*mean(mean(accel.x(ySubset,xSubset,t)));
            energy.forceY(t) = specimen.mass*mean(mean(accel.y(ySubset,xSubset,t)));

            if t == 1
                energy.impEdgeDispAvgX(t) = mean(disp.x(impEdgeY,impEdgeX,t));
                energy.impEdgeDispAvgY(t) = mean(disp.y(impEdgeY,impEdgeX,t));
                
                if EBalOpts.impEdgeXForceOnly
                    energy.impactFromAccel(t) = energy.forceX(t)*energy.impEdgeDispAvgX(t);
                else
                    energy.impactFromAccel(t) = energy.forceX(t)*energy.impEdgeDispAvgX(t)...
                    +energy.forceY(t)*energy.impEdgeDispAvgY(t);
                end
            else
                energy.impEdgeDispAvgX(t) = mean(disp.x(impEdgeY,impEdgeX,t)) - ...
                    mean(disp.x(impEdgeY,impEdgeX,t-1));
                energy.impEdgeDispAvgY(t) = mean(disp.y(impEdgeY,impEdgeX,t)) - ...
                    mean(disp.y(impEdgeY,impEdgeX,t-1));
                
                if EBalOpts.impEdgeXForceOnly
                    energy.impactFromAccel(t) = energy.impactFromAccel(t-1)+...
                        energy.forceX(t)*energy.impEdgeDispAvgX(t); 
                else
                   energy.impactFromAccel(t) = energy.impactFromAccel(t-1)+...
                        energy.forceX(t)*energy.impEdgeDispAvgX(t)+...
                        energy.forceY(t)*energy.impEdgeDispAvgY(t); 
                end
            end
        end
        
        %------------------------------------------------------------------
        % KINETIC ENERGY
        % Calculate the kinetic energy components
        % rho = m/V, m = rho*V
        energy.kineticX(t) = sum(sum(0.5*element.mass.*...
            squeeze(vel.x(ySubset,xSubset,t)).^2));
        energy.kineticY(t) = sum(sum(0.5*element.mass.*...
            squeeze(vel.y(ySubset,xSubset,t)).^2));
        energy.kinetic(t) = energy.kineticX(t)+energy.kineticY(t);
        
        %------------------------------------------------------------------
        % STRAIN ENERGY
        % Calculate the strain energy components and sum them
        energy.strainXX(t) = sum(sum(0.5*element.volume*squeeze(...
            stress.x(ySubset,xSubset,t).*strain.x(ySubset,xSubset,t))));
        energy.strainYY(t) = sum(sum(0.5*element.volume*squeeze(...
            stress.y(ySubset,xSubset,t).*strain.y(ySubset,xSubset,t))));
        energy.strainXY(t) = sum(sum(0.5*element.volume*squeeze(...
            stress.s(ySubset,xSubset,t).*strain.s(ySubset,xSubset,t))));
        energy.strain(t) = energy.strainXX(t)+energy.strainYY(t)+energy.strainXY(t);
        
    end

end

