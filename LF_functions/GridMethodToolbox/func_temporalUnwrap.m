function phases = func_temporalUnwrap(phi,threshold,method,range)
% This function temporally unwraps a sequence of 2D phases maps
% CAUTION: assumes that change in phases is small (i.e. << grid pitch)
% between successive images
%
% Inputs: 
% phi - data struct with 'x' and 'y' fields
% phi.x - 3D array containing the x phases over the 2D image, 3rd dimension
% is the time
% phi.y - 3D array containing the y phases over the 2D image, 3rd dimension
% is the time
% threshold - limit used to detect discontinuous jumps in the phases
%
% Outputs:
% phases - struct containing the following fields
% phases.x - vector with a frame by frame integer shift
% phases.y - vector with a frame by frame integer shift

if nargin < 3
    method = 'point';
elseif nargin < 4 && strcmp(method,'field')
    range.x = 1:size(phi.x,2);
    range.y = 1:size(phi.x,1);
end

numFrames = size(phi.x,3);

if strcmp(method,'field')
    phiTestX = squeeze(mean(mean(phi.x(range.y,range.x,:))));
    phiTestY = squeeze(mean(mean(phi.y(range.y,range.x,:))));       
else
    % Reference point used to detect jumps
    % auto picks a point near the centre of the image
    refPtY = round(size(phi.x,1)/2);
    refPtX = round(size(phi.x,2)/2);
    phiTestX = squeeze(phi.x(refPtY,refPtX,:));
    phiTestY = squeeze(phi.y(refPtY,refPtX,:));   
end


% Calculate diff between elements to detect jumps
stepX = diff(phiTestX);
stepY = diff(phiTestY);

phasesCountX = 0;
phasesCountY = 0;
for i = 1:numFrames-1
   if abs(stepX(i)) > threshold;
       stepSize = round(stepX(i)/(2*pi)); 
       phasesCountX = phasesCountX-stepSize;
   end

   if abs(stepY(i)) > threshold;
       stepSize = round(stepY(i)/(2*pi)); 
       phasesCountY = phasesCountY-stepSize;
   end
   phases.x(i) = phasesCountX;
   phases.y(i) = phasesCountY;
end

% Pad with a zero because the first phases map is not shifted
phases.x = [0,phases.x];
phases.y = [0,phases.y];

end

