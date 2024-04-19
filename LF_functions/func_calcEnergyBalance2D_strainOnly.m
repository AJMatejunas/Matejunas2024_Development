function energy = func_calcEnergyBalance2D_strainOnly(EBalOpts,pos,time,specimen,...
    material,element,disp,vel,accel,strain,stress)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 8/8/2017
% Calculates the kinetic, strain and impact energy components using input
% variables.
    
    % Check if we have the stress field, if not we need to get it for each
    % frame... this saves heaps of RAM though.
    calcStress = false;
    if nargin < 11
        calcStress = true;
    end

    % Remove the options from the struct to shorten code length
    sv = EBalOpts.sv;
    xSubset = EBalOpts.xSubset;
    xBound = EBalOpts.xBound;
    ySubset = EBalOpts.ySubset;
    yBound = EBalOpts.yBound;
    impEdgeX = EBalOpts.impEdgeX;
    impEdgeY = EBalOpts.impEdgeY;
    
    for t = 1:time.numFrames
        if calcStress
            stress = func_calcStressFromStrain2D_SingleFrame(strain,material.Q,t);          
        end
        
        %------------------------------------------------------------------
        % EXTERNAL WORK - Calculate from constitutive law 
        EI = func_calcImpactEnergyStressDisp(EBalOpts,pos,specimen,disp,stress,t);
        if t == 1
            energy.impactStrain(t) = sum(EI);
            clear EI      
        else
            if EBalOpts.impactEnerCutOff && (t>EBalOpts.impactEnerCutOffFrame)
                energy.impactStrain(t) = energy.impactStrain(t-1);
            else
                energy.impactStrain(t) = energy.impactStrain(t-1)+sum(EI);
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
                    energy.impactAccel(t) = energy.forceX(t)*energy.impEdgeDispAvgX(t);
                else
                    energy.impactAccel(t) = energy.forceX(t)*energy.impEdgeDispAvgX(t)...
                    +energy.forceY(t)*energy.impEdgeDispAvgY(t);
                end
            else
                if EBalOpts.impactEnerCutOff && (t>EBalOpts.impactEnerCutOffFrame)
                    energy.impactAccel(t) = energy.impactAccel(t-1);
                else
                    energy.impEdgeDispAvgX(t) = mean(disp.x(impEdgeY,impEdgeX,t)) - ...
                        mean(disp.x(impEdgeY,impEdgeX,t-1));
                    energy.impEdgeDispAvgY(t) = mean(disp.y(impEdgeY,impEdgeX,t)) - ...
                        mean(disp.y(impEdgeY,impEdgeX,t-1));

                    if EBalOpts.impEdgeXForceOnly
                        energy.impactAccel(t) = energy.impactAccel(t-1)+...
                            energy.forceX(t)*energy.impEdgeDispAvgX(t); 
                    else
                       energy.impactAccel(t) = energy.impactAccel(t-1)+...
                            energy.forceX(t)*energy.impEdgeDispAvgX(t)+...
                            energy.forceY(t)*energy.impEdgeDispAvgY(t); 
                    end
                    
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
        if calcStress
            energy.strainXX(t) = sum(sum(0.5*element.volume*squeeze(...
                stress.x(ySubset,xSubset).*strain.x(ySubset,xSubset,t))));
            energy.strainYY(t) = sum(sum(0.5*element.volume*squeeze(...
                stress.y(ySubset,xSubset).*strain.y(ySubset,xSubset,t))));
            energy.strainXY(t) = sum(sum(0.5*element.volume*squeeze(...
                stress.s(ySubset,xSubset).*strain.s(ySubset,xSubset,t))));
        else
            energy.strainXX(t) = sum(sum(0.5*element.volume*squeeze(...
                stress.x(ySubset,xSubset,t).*strain.x(ySubset,xSubset,t))));
            energy.strainYY(t) = sum(sum(0.5*element.volume*squeeze(...
                stress.y(ySubset,xSubset,t).*strain.y(ySubset,xSubset,t))));
            energy.strainXY(t) = sum(sum(0.5*element.volume*squeeze(...
                stress.s(ySubset,xSubset,t).*strain.s(ySubset,xSubset,t))));
        end
        energy.strain(t) = energy.strainXX(t)+energy.strainYY(t)+2*energy.strainXY(t);
        
    end
    
    %----------------------------------------------------------------------
    % Calculate the energy balance 
    energy.kineticPlusStrain = energy.kinetic + energy.strain;
    
    % Impact Energy From Strain
    energy.balStrain = energy.impactStrain - energy.kinetic - energy.strain;
    energy.balStrainErrPc = max(abs(energy.balStrain(EBalOpts.maxErrRange)))...
        /max (abs(energy.impactStrain(EBalOpts.maxErrRange)))*100;

    % Impact Energy From Accel
    if EBalOpts.calcExtWorkFromAccel
        energy.balAccel = energy.impactAccel - energy.kinetic - energy.strain;
        energy.balAccelErrPc = max(abs(energy.balAccel(EBalOpts.maxErrRange)))...
            /max (abs(energy.impactAccel(EBalOpts.maxErrRange)))*100;
    end
end

