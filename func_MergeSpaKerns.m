function [MergedVec,MergedVar]=func_MergeSpaKerns(SmallVec,...
    LargeVec,SmallVar,LargeVar)

%This function is written to combine identification results for different
    %smoothing sweeps

% Author: Andrew Matejunas

%Date Created: 2023/04/07


% Function input arguments
    %SmallVec- Vector of smaller spatial smoothing kernels
    %LargeVec- Vector of larger spatial smoothing kernels
    %SmallVar- Matrix containing variable to be combined for smaller
                    %kernels [TemporalSmoothing,SpatialSmoothing]
    %LargeVar- Matrix containing variable to be combined for larger kernels
                    %[TemporalSmoothing,SpatialSmoothing]


% Function output arguments
    %MergedVec- vector of merged smoothing kernels
    %MergedVar- merged matrix of the variable being combined
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                   
                    
MergedVec=[SmallVec,LargeVec];
MergedVar=[SmallVar,LargeVar];
end