function  disp = func_fixNaNsInDisp(extrapOpts,disp)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 20/6/2017
%
% Finds NaNs in the displacement fields and replace them with local
% averages

% Find the location of all the NaNs in the arrays
[y,x,t] = ind2sub(size(disp),find(isnan(disp)));
numNans = length(x);

for i = 1:numNans
   % Set the ranges for the mean replacement value for the NaN
   xRange = x(1)-extrapOpts.fixNaNKernal:x(1)+extrapOpts.fixNaNKernal;
   % Make sure we are not out of bounds
   if min(xRange) < 1
       xRange = 1:x(1)+extrapOpts.fixNaNKernal;
   elseif max(xRange) > size(disp,2)
       xRange = x(1)-extrapOpts.fixNaNKernal:size(disp,2); 
   end
   yRange = y(1)-extrapOpts.fixNaNKernal:y(1)+extrapOpts.fixNaNKernal;
   if min(yRange) < 1
       yRange = 1:y(1)+extrapOpts.fixNaNKernal;
   elseif max(yRange) > size(disp,1)
       yRange = y(1)-extrapOpts.fixNaNKernal:size(disp,1); 
   end
   
   checkMean = squeeze(nanmean(nanmean(disp(yRange,xRange,t(i)))));
   if isnan(checkMean)
       fprintf('\t\tWARNING: NaN replacement at [%i,%i,%i] is still NaN, widen kernal.\n',x(i),y(i),t(i))
   end
   disp(y(i),x(i),t(i)) = checkMean;  
end

end

