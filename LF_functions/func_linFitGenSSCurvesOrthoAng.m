function identStiffGenSS = func_linFitGenSSCurvesOrthoAng(genSSCurveOpts,...
    angledSlice1,angledSlice2,genSSCurves)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 29/1/2019
% Date Edited: 29/1/2019
%
% Linearly fits the generalised stress-strain curve formulas to identify
% the stiffness components for orthotropic elasticity: Q11, Q22, Q12, Q66

    if angledSlice1.xMaxInd > 0
        % Calculate indices to average over the middle slices
        genSSCurveOpts.slice1AvgQVsLRange = round(genSSCurveOpts.avgQVsLRangePc(1)...
            *angledSlice1.xMaxInd):round(genSSCurveOpts.avgQVsLRangePc(2)*angledSlice1.xMaxInd);
        % Fit the stress strain curves to obtain Q22, Q12 and Q66 using 
        % equations 14.1 to 14.3
        % Eq14.1 Q22
        fprintf('\tLinearly fitting slice 1 for Q22.\n')
        [identStiffGenSS.Q22VsL_Eq141,~] = ...
            func_identStiffLinFitStressStrainCurve(genSSCurveOpts,...
            genSSCurves.Q22_accelAvg_Eq141,genSSCurves.Q22_strainAvg_Eq141);
        identStiffGenSS.Q22AvgOverL_Eq141 = nanmean(...
            identStiffGenSS.Q22VsL_Eq141(genSSCurveOpts.slice1AvgQVsLRange));

        % Eq14.2 Q12
        fprintf('\tLinearly fitting slice 1 for Q12.\n')
        [identStiffGenSS.Q12VsL_Eq142,~] = ...
            func_identStiffLinFitStressStrainCurve(genSSCurveOpts,...
            genSSCurves.Q12_accelAvg_Eq142,genSSCurves.Q12_strainAvg_Eq142);
        identStiffGenSS.Q12AvgOverL_Eq142 = nanmean(...
            identStiffGenSS.Q12VsL_Eq142(genSSCurveOpts.slice1AvgQVsLRange));

        % Eq14.3 Q66
        fprintf('\tLinearly fitting slice 1 for Q66.\n')
        [identStiffGenSS.Q66VsL_Eq143,~] = ...
            func_identStiffLinFitStressStrainCurve(genSSCurveOpts,...
            genSSCurves.Q66_accelAvg_Eq143,genSSCurves.Q66_strainAvg_Eq143);
        identStiffGenSS.Q66AvgOverL_Eq143 = nanmean(...
            identStiffGenSS.Q66VsL_Eq143(genSSCurveOpts.slice1AvgQVsLRange));
    end

    if angledSlice2.xMaxInd > 0
        % Calculate indices to average over the middle slices
        genSSCurveOpts.slice2AvgQVsLRange = round(genSSCurveOpts.avgQVsLRangePc(1)...
            *angledSlice2.xMaxInd):round(genSSCurveOpts.avgQVsLRangePc(2)*angledSlice2.xMaxInd);

        % Fit the stress strain curves to obtain Q11, Q12 and Q66 using 
        % equations 18.1 to 18.3
        % Eq18.1 Q11
        fprintf('\tLinearly fitting slice 2 for Q11.\n')
        [identStiffGenSS.Q11VsL_Eq181,~] = ...
            func_identStiffLinFitStressStrainCurve(genSSCurveOpts,...
            genSSCurves.Q11_accelAvg_Eq181,genSSCurves.Q11_strainAvg_Eq181);
        identStiffGenSS.Q11AvgOverL_Eq181 = nanmean(...
            identStiffGenSS.Q11VsL_Eq181(genSSCurveOpts.slice2AvgQVsLRange));

        % Eq18.2 Q12
        fprintf('\tLinearly fitting slice 2 for Q12.\n')
        [identStiffGenSS.Q12VsL_Eq182,~] = ...
            func_identStiffLinFitStressStrainCurve(genSSCurveOpts,...
            genSSCurves.Q12_accelAvg_Eq182,genSSCurves.Q12_strainAvg_Eq182);
        identStiffGenSS.Q12AvgOverL_Eq182 = nanmean(...
            identStiffGenSS.Q12VsL_Eq182(genSSCurveOpts.slice2AvgQVsLRange));

        % Eq18.3 Q66
        fprintf('\tLinearly fitting slice 2 for Q66.\n')
        [identStiffGenSS.Q66VsL_Eq183,~] = ...
            func_identStiffLinFitStressStrainCurve(genSSCurveOpts,...
            genSSCurves.Q66_accelAvg_Eq183,genSSCurves.Q66_strainAvg_Eq183);
        identStiffGenSS.Q66AvgOverL_Eq183 = nanmean(...
            identStiffGenSS.Q66VsL_Eq183(genSSCurveOpts.slice2AvgQVsLRange));
    end
end

