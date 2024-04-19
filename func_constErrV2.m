function [Errors]=func_constErrV2(ABQstress,...
            ConstStress)

% This function is written to compute the errors between the stress output
    % by abaqus for an IBII experiment and the stress calculated through
    % the constitutive model

%inputs
    %ABQstress- full structure of ABAQUS stresses
    %ConstStress- full structure of stresses calculated with the
        %constitutive model
%outputs
    %stressXXerror- error in the normal stress in the X direction
    %stressXYerror- error in the shear stresss
    %stressYYerror- error in normal stress in Ydirection
    %XXdiff,YYdiff,Sheardiff- differences between constitutive stress and
        %Abaqus stress
    
%% Reorder the inputs into an easier form for calculation

ABQX=ABQstress.x;
ABQXY=ABQstress.s;
ABQY=ABQstress.y

MaxABQx=max(abs(ABQX),[],'all');
MaxABQxy=max(abs(ABQXY),[],'all');
MaxABQy=max(abs(ABQY),[],'all');

ConstX=ConstStress.xx;
ConstY=ConstStress.yy;
ConstXY=ConstStress.xy;

%%Calculate Differences
XXdiff=ConstX-ABQX;
YYdiff=ConstY-ABQY;
Sheardiff=ConstXY-ABQXY;

%% calculate the errors

stressXXerr=(XXdiff)./ABQX*100;
stressXXerr(ABQX==0)=0;
stressYYerr=(YYdiff-ABQY)./ABQY*100;
stressYYerr(ABQY==0)=0;
stressXYerr=(Sheardiff)./ABQXY*100;
stressXYerr(ABQXY==0)=0;

%% calculate Normalized Errors
XXnormErr=abs(XXdiff./MaxABQx)*100;
XYnormErr=abs(XYdiff./MaxABQxy)*100;
YYnormErr=abs(YYdiff./MaxABQy)*100;
end
