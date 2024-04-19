function [outputArg1,outputArg2] = func_ProcessDICDisp(specimen,...
    material,)
%This code is written to process grid method displacements from MatchID.
%Note this version of the code does not extract strains. Strains are
%calculated later in the same manner as the grid method. 




globalOpts.ExtractDICStrain=false;


outputArg1 = inputArg1;
outputArg2 = inputArg2;
end