function [avgVar] = func_avgFFVarOverWidth(ffVar)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 27/3/2017
%
% Averages full field data over the width, assumes the first dimension of
% the array represents the transverse or 'Y' direction and averages over
% this dimension.

avgVar = squeeze(nanmean(ffVar,1));

end

