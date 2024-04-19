function [extrapVec,r] = func_quick1DLinearExtrap(xToFit,dispToFit,...
    xEdge,dispEdge,xToExtrap)
% Quick 1D linear extrapolation function, used for extrapolating missing
% edges data on displacement fields from the grid method or DIC.
%
% Author: Jared Van-Blitterswyk, Lloyd Fletcher
% Image-Based Mechanics Group (IBMG), University of Southampton
% Date Created: 28/7/2020
% Date Edited: 28/7/2020

    x = [ones(length(xToFit),1) xToFit'];        
    r = x\dispToFit';
    diff = dispEdge - (r(2)*xEdge + r(1));
    r(1) = r(1) + diff;
    extrapVec = r(2)*xToExtrap + r(1);
end

