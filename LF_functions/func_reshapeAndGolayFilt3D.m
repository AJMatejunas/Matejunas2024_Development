function filteredVar = func_reshapeAndGolayFilt3D(varToFilt,polyOrder,frameSize)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 29/8/2017
% Date Edited: 30/7/2020 - updated to not require the smoothingOpts structure 

%Moved to old functions folder 2023/03/04

% Flatten spatial dimensions and filter with savitsky-golay filter.
% Required for forward compatibility as 2017 version of matlab returns a
% flat array when given a 3D matrix.

[sy,sx,st] = size(varToFilt);

% Flatten the 3D array into 2D
flatVar = reshape(varToFilt,[sy*sx,st]);

% Filter the flat array
filteredVar = sgolayfilt(flatVar,polyOrder,frameSize,[],2);

% Reshape back to 3D and return this parameter
filteredVar = reshape(filteredVar,[sy,sx,st]);

end

