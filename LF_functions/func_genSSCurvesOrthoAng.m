function genSSCurves = func_genSSCurvesOrthoAng(material,angledSlice1,angledSlice2,...
    strainSlice1,strainSlice2,accelSlice1,accelSlice2,debug)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 29/1/2019
% Date Modified: 29/1/2019 
%
% Calculates stress metrics for the generalised stress-strain curves in
% orthotropic elasticity. These equations give a formulation that allows
% various weighted averages to be plotted against each other to identify
% the in plane stiffness components Q11, Q12, Q22 and Q66. 
%
% Full details can be found in the following paper:
% Pierron, F. and Fletcher, L., "Image-based generalized stress-strain 
% curves for inertial high strain rate tests on isotropic and orthotropic
% materials" 
%
% The function takes input data structs describing the kinematic fields in
% specimen and outputs a data structure containing the values of each
% equation set allowing the linear relationship for each Q component to be 
% plotted and fitted to obtain the relevant stiffness component.
    
    if nargin < 8
        debug = false;
    end

    % Section to plot the 'test curve', pick the middle of the sample
    tc1 = round(angledSlice1.xMaxInd/2);
    tc2 = round(angledSlice2.xMaxInd/2);
    
    % Equations 14 - 1,2,3
    if angledSlice1.xMaxInd > 0
        % Eq 14.1
        genSSCurves.Q22_strainAvg_Eq141 = strainSlice1.avg11.*strainSlice1.avg22_x1 - ...
            strainSlice1.avg22.*strainSlice1.avg11_x1;
        genSSCurves.Q22_accelAvg_Eq141 = material.rho*angledSlice1.surfArea./angledSlice1.length.*...
            ((accelSlice1.surfAvg11_x2 - accelSlice1.surfAvg22_x1).*strainSlice1.avg11 +...
            accelSlice1.surfAvg22.*strainSlice1.avg11_x1);

        % Eq 14.2
        genSSCurves.Q12_strainAvg_Eq142 = strainSlice1.avg22.*strainSlice1.avg11_x1 - ...
            strainSlice1.avg11.*strainSlice1.avg22_x1;
        genSSCurves.Q12_accelAvg_Eq142 = material.rho*angledSlice1.surfArea./angledSlice1.length.*...
            ((accelSlice1.surfAvg11_x2 - accelSlice1.surfAvg22_x1).*strainSlice1.avg22 +...
            accelSlice1.surfAvg22.*strainSlice1.avg22_x1);

        % Eq 14.3
        genSSCurves.Q66_strainAvg_Eq143 = strainSlice1.avg12;
        genSSCurves.Q66_accelAvg_Eq143 = -material.rho*angledSlice1.surfArea./angledSlice1.length.*...
            (accelSlice1.surfAvg11);
        
        if debug
            eqSetStr = 'Eq. 14';
            paramStr = 'Q_{22}';
            tStr = ['Stress-Strain Curve, ',paramStr,': ',eqSetStr];
            X = squeeze(genSSCurves.Q22_strainAvg_Eq141(tc1,:))';
            Y = squeeze(genSSCurves.Q22_accelAvg_Eq141(tc1,:))';
            hf = func_plotTestSSCurve(paramStr,tStr,X,Y);
            
            paramStr = 'Q_{12}';
            tStr = ['Stress-Strain Curve, ',paramStr,': ',eqSetStr];
            X = squeeze(genSSCurves.Q12_strainAvg_Eq142(tc1,:))';
            Y = squeeze(genSSCurves.Q12_accelAvg_Eq142(tc1,:))';
            hf = func_plotTestSSCurve(paramStr,tStr,X,Y);
            saveFile = [savePath,'\','Eq_14_2'];
            print(hf,saveFile,plotProps.format,plotProps.saveRes)
            saveas(hf,saveFile,'fig')
        
            paramStr = 'Q_{66}';
            tStr = ['Stress-Strain Curve, ',paramStr,': ',eqSetStr];
            X = squeeze(genSSCurves.Q66_strainAvg_Eq143(tc1,:))';
            Y = squeeze(genSSCurves.Q66_accelAvg_Eq143(tc1,:))';
            hf = func_plotTestSSCurve(paramStr,tStr,X,Y);
        end
    end

    % Equation 15
    if angledSlice1.xMaxInd > 0
        genSSCurves.Q22_strainAvg_Eq15 = strainSlice1.avg22;
        genSSCurves.Q22_accelAvg_Eq15 = -material.rho*angledSlice1.surfArea./...
                angledSlice1.length.*accelSlice1.surfAvg22;
        
        if debug
            eqSetStr = 'Eq. 15';
            % Plot a test stress strain curve
            paramStr = 'Q_{22}';
            tStr = ['Stress-Strain Curve, ',paramStr,': ',eqSetStr];
            X = squeeze(strainSlice1.avg22(tc1,:))';
            Y = squeeze(stressSlice1.avg22(tc1,:))';
            hf = func_plotTestSSCurve(paramStr,tStr,X,Y);
        end
    end

    % Equations 18 - 1,2,3
    if angledSlice2.xMaxInd > 0
        % Eq 18.1
        genSSCurves.Q11_strainAvg_Eq181 = strainSlice2.avg22.*strainSlice2.avg11_x2 - ...
           strainSlice2.avg11.*strainSlice2.avg22_x2;
        genSSCurves.Q11_accelAvg_Eq181 = material.rho*angledSlice2.surfArea./angledSlice2.length.*...
            ((accelSlice2.surfAvg11_x2 - accelSlice2.surfAvg22_x1).*strainSlice2.avg22 -...
            accelSlice2.surfAvg11.*strainSlice2.avg22_x2);

        % Eq 18.2
        genSSCurves.Q12_strainAvg_Eq182 = strainSlice2.avg11.*strainSlice2.avg22_x2 - ...
           strainSlice2.avg22.*strainSlice2.avg11_x2;
        genSSCurves.Q12_accelAvg_Eq182 = material.rho*angledSlice2.surfArea./angledSlice2.length.*...
            ((accelSlice2.surfAvg11_x2 - accelSlice2.surfAvg22_x1).*strainSlice2.avg11 -...
            accelSlice2.surfAvg11.*strainSlice2.avg11_x2);

        % Eq 18.3
        genSSCurves.Q66_strainAvg_Eq183 = strainSlice2.avg12;
        genSSCurves.Q66_accelAvg_Eq183 = material.rho*angledSlice2.surfArea./angledSlice2.length.*...
            (accelSlice2.surfAvg22);
        
        if debug
            eqSetStr = 'Eq. 18';
            paramStr = 'Q_{11}';
            tStr = ['Stress-Strain Curve, ',paramStr,': ',eqSetStr];
            X = squeeze(genSSCurves.Q11_strainAvg_Eq181(tc2,:))';
            Y = squeeze(genSSCurves.Q11_accelAvg_Eq181(tc2,:))';
            hf = func_plotTestSSCurve(paramStr,tStr,X,Y);

            paramStr = 'Q_{12}';
            tStr = ['Stress-Strain Curve, ',paramStr,': ',eqSetStr];
            X = squeeze(genSSCurves.Q12_strainAvg_Eq182(tc2,:))';
            Y = squeeze(genSSCurves.Q12_accelAvg_Eq182(tc2,:))';
            hf = func_plotTestSSCurve(paramStr,tStr,X,Y);

            paramStr = 'Q_{66}';
            tStr = ['Stress-Strain Curve, ',paramStr,': ',eqSetStr];
            X = squeeze(genSSCurves.Q66_strainAvg_Eq183(tc1,:))';
            Y = squeeze(genSSCurves.Q66_accelAvg_Eq183(tc1,:))';
            hf = func_plotTestSSCurve(paramStr,tStr,X,Y);
        end
    end
   
end



