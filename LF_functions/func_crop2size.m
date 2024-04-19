function [Cfield] = func_crop2size(field,Rx,Ry)
% This function is written to crop an extrapolated kinematic field back to
% its original size if the extapolation window is larger than the ROI

%inputs
    %field- extrapolated kinematic field

%outputs
    %Cfield- field cropped back to original size


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cfield=field(Ry,Rx,:);

end