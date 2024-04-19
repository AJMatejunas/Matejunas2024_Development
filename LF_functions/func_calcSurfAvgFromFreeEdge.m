function [ surfAvg ] = func_calcSurfAvgFromFreeEdge( inputField )
% Author: L. Fletcher
% PhotoDyn Research Group
% Date Created: 3/12/2018 
% Date Edited: 3/12/2018
%
% Calcultes a surface average of a variable from the free edge to the given
% row, similar to the stress-gauge equation.

% Get the size of the input matrix
[sy,sx,st] = size(inputField);

surfAvg = zeros([sx,st]);
for tt = 1:st
    for xx = 1:sx
        tempField = inputField(:,1:xx,tt);
        surfAvg(xx,tt) = mean(tempField(:));
    end
end


end

