function [RMSE] = func_calcFFRMSE(RefField,SimField)
%This function is written to calculate the full field RMSE error given a
    %reference data (Generally FE kinematic fields interpolated to Grid
    %Method) and kinematic fields obtained by the grid method

%Author: Andrew Matejunas

%Date: 2023/03/08

%Version History/Chage Log
 
%Function Input arguments
    %RefField- Array containing the reference values
    %SimField- Array containig values obtained through simulations

%Function outputs:
    %RMSE- Root mean square error in the same units as the field defined as
        %sqrt(1/numberOfPoints*(Sim-Ref)^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Determine number of elements
TotNum=numel(RefField);

%% Calculate Difference in field value at each point
Error=(SimField-RefField);

%% Calculate RMSE
RMSE=sqrt((1/TotNum)*sum(Error.^2,'all'));

end