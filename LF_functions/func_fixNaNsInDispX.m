function  disp = func_fixNaNsInDispX(extrapOpts,disp)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 03/5/2018
%
% Finds NaNs in the displacement fields and replace them with local
% averages
nKernel = extrapOpts.fixNaNKernel; % total width of kernel 
% Find the location of all the NaNs in the arrays
[y,x,t] = ind2sub(size(disp),find(isnan(disp)));
numNans = length(x);

for i = 1:numNans
    
   xRange = x(i)-ceil(0.5*nKernel):x(i)+ceil(0.5*nKernel);
    if min(xRange) < 1
        xRange = 1:x(i)+ceil(0.5*nKernel);
    elseif max(xRange) > size(disp,2)
       xRange = x(i)-ceil(0.5*nKernel):size(disp,2); 
    end
    
   yRange = y(i)-ceil(0.5*nKernel):y(i)+ceil(0.5*nKernel);
   if min(yRange) < 1
       yRange = 1:y(i)+ceil(0.5*nKernel);
   elseif max(yRange) > size(disp,1)
       yRange = y(i)-ceil(0.5*nKernel):size(disp,1); 
   end
   
   disp(y(i),x(i),t(i)) = squeeze(nanmean(nanmean(disp(yRange,xRange,t(i)))));  
end

end

