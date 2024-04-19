function [xxSGt,ShearSGt,ConstStressXt,ConstStressSt,CondOpts] =...
    func_thresholdKinFields(xxSG,ShearSG,ConstStressX,ConstStressS,...
    strainXX,strainS,CondOpts)
%This function is written to threshold the stress fields based on the
    %measured resolution in the strain fields and stress gauge stresses.
    %the noise floor is measured from undeformed images, and this code is
    %used to remove data below a selected threshold value from the input
    %into the cost function.

%Note this version of the code exludes all data that is below the noise
    %floor from the calulation of phi. It may be discovered later that
    %either an initial point of 0 strain 0 stress, or some data might be
    %needed during relaxation

    
%Author: Andrew Matejunas

%Date Created: 2023/01/11

%Version History/Change log:

%Function Input Arguments:
    %xxSG- Matrix of normal stress gauge stresses in the x-direction
        %[x-locations,time points] (Pa)
    %ShearSG- Matrix of stress gauge in-plane shear stresses 
        %[x-locations,time points] (Pa)
    %ConsStressX- Matrix of average normal stresses in the x-direction
        %calculated from the constitutive model [x-locations,time points]
        %(Pa)
    %ConsStressS- Matrix of average in-plane shear stresses calculated from
        %the constitutive model [x-locations,time points] (Pa)
    %strainXX- Matrix of average normal strains in the x-direction
        % [x-locations,time points] used to threshold constitutive model
        % stresses are above the noise floor
    %strainS- Matrix of average in-plane shear strains 
        % [x-locations,time points] used to threshold constitutive model
        % stresses are above the noise floor
    %CondOpts- structure containing the data processing options. The 
        %critical fields for this particulat function include
            %SGxNF- noise floor in the stress gauge stresses in the
                %x-direction
            %SGsNF- noise floor in the shear stress gauge stresses
            %StrainxNF- noise floor in the normal strains in the
                %x-direction
            %StrainsNF- noise floor in the in-plane shear strains
            %ThresFact- multiplaction factor to determine how far above the
                %noise floor to censor the data

%Function output arguments:
    %xxSG- Matrix of thresholded normal stress gauge stresses in the 
        %x-direction [x-locations,time points] (Pa)
    %ShearSGt- Matrix of thresholded stress gauge in-plane shear stresses 
        %[x-locations,time points] (Pa)    
    %ConsStressXt- Matrix of thresholded average normal stresses in the
        %x-direction calculated from the constitutive model 
        %[x-locations,time points] (Pa)
    %ConsStressSt- Matrix of thresholded average in-plane shear stresses 
        %calculated from the constitutive model [x-locations,time points]
        %(Pa)
    %CondOpts- structure containing the data processing options. The 
        %critical fields for this particulat function include
            %SGxNF- noise floor in the normal stress gauge stress 
            %SGsNF- noise floor in the shear stress gauge stress
            %StrainxNF- noise floor in the normal strains in the
                %x-direction
            %StrainsNF- noise floor in the in-plane shear strains
            %ThresFact- multiplaction factor to determine how far above the
                %noise floor to censor the data
            %Thresh- structure containing the actual thresholding limits



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Determine thresholding Parmaters
Thresh.SGx=CondOpts.ThreshFact*CondOpts.SGxNF;
Thresh.SGs=CondOpts.ThreshFact*CondOpts.SGsNF;
Thresh.Strainx=CondOpts.ThreshFact*CondOpts.StrainxNF;
Thresh.Strains=CondOpts.ThreshFact*CondOpts.StrainxNF;

%% Vecorize kinematic fields
vecSize=numel(xxSG);
SGvecX=reshape(xxSG,[vecSize,1]);
SGvecS=reshape(ShearSG,[vecSize,1]);
ConstVecX=reshape(ConstStressX,[vecSize,1]);
ConstVecS=reshape(ConstStressS,[vecSize,1]);
StrainVecX=reshape(strainXX,[vecSize,1]);
StrainVecS=reshape(strainS,[vecSize,1]);

%% Threshold data
%Only considering data where the average strain and 
xxSGt=SGvecX(abs(SGvecX)>=Thresh.SGx & abs(StrainVecX)>=Thresh.Strainx);
ConstStressXt=ConstVecX(abs(SGvecX)>=...
    Thresh.SGx & abs(StrainVecX)>=Thresh.Strainx);
ShearSGt=SGvecS(abs(SGvecS)>=Thresh.SGs & abs(StrainVecS)>=Thresh.StrainS);
ConstStressSt=ConstVecS(abs(SGvecS)>=...
    Thresh.SGs & abs(StrainVecS)>=Thresh.Strains);
%% Store outputs
CondOpts.Thresh=Thresh;
end