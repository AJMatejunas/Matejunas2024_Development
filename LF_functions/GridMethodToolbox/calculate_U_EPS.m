function [ UX, UY, EPSXX, EPSYY, EPSXY, WXY ] = calculate_U_EPS( p, PHI1_X, PHI1_Y, PHI2_X, PHI2_Y, procedure, maxiter, suppressOutput)
% calculate displacement and strain components from phase maps
% 1) input:
% p: grid pitch 
% PHI1_X, PHI1_Y, PHI2_X, PHI2_Y: phase maps in non-deformed (1) and deformed (2) grid images
% procedure: 1 = approximate, 2=PHI2 backdeformed in PHI1 axis, 2=iterative (until relative improvement below "prec" value)
% maxiter: for procedure 2, maximum number of iterations until relative improvement below "stopcrit" value (optional parameter, default: maxiter=1)
% 2) output:
% UX, UY: displacement maps
% EPSXX, EPSYY, EPXY: strain maps
% WXY: rotation map


if nargin == 6
    maxiter=50; 
end

if nargin < 8
    suppressOutput = false;
end

UX = -( PHI2_X-PHI1_X ) * p/(2*pi);
UY = -( PHI2_Y-PHI1_Y ) * p/(2*pi);

if (procedure==2)
    [Y,X]=ndgrid(1:size(PHI1_X,1),1:size(PHI1_X,2));
    stopcrit=5e-4; % smaller values useless because of interpolation bias
    for n=1:maxiter;
        interpPHI2_X= interp2(X,Y,PHI2_X,X+UX,Y+UY,'cubic');
        interpPHI2_Y= interp2(X,Y,PHI2_Y,X+UX,Y+UY,'cubic');
        newUX = -( interpPHI2_X-PHI1_X ) * p/(2*pi);
        newUY = -( interpPHI2_Y-PHI1_Y ) * p/(2*pi);
        deltaX=(newUX(:)-UX(:)); deltaY=(newUY(:)-UY(:));
        UX=newUX;
        UY=newUY;
        if ((norm(deltaX(isfinite(deltaX)))<stopcrit*numel(deltaX))&&(norm(deltaY(isfinite(deltaY)))<stopcrit*numel(deltaY))) 
            break; 
        end
    end
    if ~suppressOutput
        disp(['procedure 2: ',num2str(n),' iteration(s)  -- delta x:',num2str(norm(deltaX(isfinite(deltaX)))/numel(deltaX)),'  delta y:',num2str(norm(deltaX(isfinite(deltaY)))/numel(deltaY))]);
    end
end

[GxUX,GyUX]=gradient(UX);
[GxUY,GyUY]=gradient(UY);

EPSXX=GxUX;
EPSYY=GyUY;
EPSXY=(GyUX+GxUY)/2;        
WXY=(GyUX-GxUY)/2; 

end

