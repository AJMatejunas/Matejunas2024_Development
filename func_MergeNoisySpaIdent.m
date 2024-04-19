function MergedIdent = func_MergeNoisySpaIdent(SmallIdent,LargeIdent)
%This function is written to combine the full noisy identification arrays
    %for two different spatial smoothing sweeps.

    %Author: Andrew Matejunas
    %Date Created: 2023/04/07

% Function input arguments
    %SmallIdent- Array containing identified consitutive parameter for
                    %smaller smoothing kernels
                    %[TempSmooth,SpatialSmooth,NoiseIteration]
    %LargeIdent- Array containing identified consitutive parameter for
                    %larger smoothing kernels
                    %[TempSmooth,SpatialSmooth,NoiseIteration] 

% Output Arguments
    %MergedIdent- Array containing merged identified constitutive parameter
                    %array

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MergedIdent=cat(2,SmallIdent,LargeIdent);
end