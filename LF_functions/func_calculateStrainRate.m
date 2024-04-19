function [strainRate,strain] = func_calculateStrainRate(strain,time,smoothingOpts,returnAllComps)
% Author: Jared Van Blitterswyk
% PhotoDyn Group, University of Southampton
% Date: 23/03/2017
% Edited By: Lloyd Fletcher
% Date Edited: 14/2/2018 - add calculation for y and s components
% This function takes in previously smoothed strain maps and computes
% returns full-field strain maps in the struct: 'strainRate'

    % If unspecified only calculate the 'x' direction strain rate
    if nargin < 4
        returnAllComps = false;
    end

    if returnAllComps % Calculate the strain rate for x,y,s
        % perform full-field temporal smoothing using sliding window polynomial
        if smoothingOpts.FFTempSmooth
            strain.xtSmooth = func_reshapeAndGolayFilt3D(strain.x,smoothingOpts);
            strain.ytSmooth = func_reshapeAndGolayFilt3D(strain.y,smoothingOpts);
            strain.stSmooth = func_reshapeAndGolayFilt3D(strain.s,smoothingOpts);
        end

        % compute strain rate
        if smoothingOpts.FFTempSmooth
            [~,~,strainRate.x] = gradient(strain.xtSmooth,1,1,time.step);
            [~,~,strainRate.y] = gradient(strain.ytSmooth,1,1,time.step);
            [~,~,strainRate.s] = gradient(strain.stSmooth,1,1,time.step);
        else
            [~,~,strainRate.x] = gradient(strain.x,1,1,time.step);
            [~,~,strainRate.y] = gradient(strain.y,1,1,time.step);
            [~,~,strainRate.s] = gradient(strain.s,1,1,time.step);
        end
    else % Only calculate the x strain rate
        % perform full-field temporal smoothing using sliding window polynomial
        if smoothingOpts.FFTempSmooth
            strain.xtSmooth = func_reshapeAndGolayFilt3D(strain.x,smoothingOpts);
        end

        % compute strain rate
        if smoothingOpts.FFTempSmooth
            [~,~,strainRate.x] = gradient(strain.xtSmooth,1,1,time.step);
        else
            [~,~,strainRate.x] = gradient(strain.x,1,1,time.step);
        end
    end

end