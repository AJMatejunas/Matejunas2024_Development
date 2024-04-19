function [extrapVec,r] = func_quick1DQuadraticExtrap(xToFit,dispToFit,...
    xToExtrap)
% Quick 1D quadratic extrapolation function, used for extrapolating missing
% edges data on displacement fields from the grid method or DIC.
%
% Author: Jared Van-Blitterswyk, Lloyd Fletcher
% Image-Based Mechanics Group (IBMG), University of Southampton
% Date Created: 28/7/2020
% Date Edited: 28/7/2020

    x = [ones(length(xToFit),1) xToFit' xToFit'.^2];
    r = x\dispToFit';
    extrapVec = r(3)*xToExtrap.^2 + r(2)*xToExtrap + r(1);
end

