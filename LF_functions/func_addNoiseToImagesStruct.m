function noisyImages = func_addNoiseToImagesStruct(imageStack,imageNoise)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 2/10/2017
%
% Adds flat (grey level independent) gaussian noise to an image stack

% If there are any NaNs in the images we must find them and remove them
if sum(sum(sum(isnan(imageStack)))) > 0
    % Set the NaNs to be grey
    fprintf('WARNING: nans in image stack!')
    imageStack(isnan(imageStack)) = (2^imageNoise.bits) / 2;    
end

% Shuffle the random generator and use a truly random 'twister' generator
rng('shuffle','twister')
% Create a 3D array of gaussian noise
noise = randn(size(imageStack));
% Add the noise to the image stack
noisyImages = imageStack + imageNoise.pcNoise/100.*2^imageNoise.bits.*noise;

if imageNoise.convToUInt16
    noisyImages = uint16(noisyImages);
end
noisyImages = double(noisyImages);

end