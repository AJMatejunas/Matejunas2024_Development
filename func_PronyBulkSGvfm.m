function [SG_param,CostFunc]=func_PronyBulkSGvfm(strain,time,SG,...
    exactProps,guessvec,ubvec,lbvec,SolveOpts,minOpts,CondOpts)

%This script is written to perform a nonlinear identification on a prony
%series viscoelastic model using the stress gage formulation of the
%nonlinear virtual fields method. The code was adapted from a combined cost
%function that simultaneously solved for K and G

%Author: Andrew Matejunas
%Created:
%Validated:

%Change log:
    %2023/01/17: Adapted from func_PronySGvfmV5. Removed shear and time
                    %constant from the cost function.
                %Cost function used is now func_ViscoBulkCostSG
   %2023/01/18: Removed ability to assume a constant Poisson's ratio or
                    %deal with unknown long term moduli
         
      
% Description
    %This script identifies material parameters by performing a nonlinear
    %minimization of the cost function:
        %phi=[(Stress calculated w/ stress gage) - (stress calculated with
        %the constitutive model)]^2
%Input parameters
    %strain- strain, either measured or output by the FE model with fields
        %x- normal strain in the x direction
        %y- normal strain in the y direction
        %s- shear strains in xy plane
    %SG- stress calculated with the stress gage equation
    %exactProps- structure containing exact material properties. Contents
        %will vary depending on the assumptions made in the formulation of
        %the problem. Common parameters are
            %Kinf- Long term bulk modulus (commonly assumed to be the
                %quasi-static bulk modulus)
            %Ginf-long term shear modulus
            %nu- QS Poisson's ratio
            %Gi- Vector of identified shear moduli
            %tau- Vector of time constants associated with the shear moduli

    %guessvec- Vector of intial guesses for bulk moduli
        %[K1,K2,...,Kn]
        
           
     %ubvec- Vector of upper bound on the bulk moduil [Ki]
     %lbvec- Vector of lower bounds on the bulk moduli [Ki]
      
     %SolveOpts- solver options for the fmincon function
        %Options inculde
            %constnu- Assumes Poisson's ratio is not dependent on strain
                %rate and remains known. Possible values are 1 and
                %0. Use 1 if nu does not depend on strain rate 
            %KGsame- Assumes the same rate dependence for K and G. Possile 
                %values are true and false. Use if ki=gi for all time 
                %constants.
            %identEinf- Determines whether E inf is assumed or a variable
                %to solve for. Possible values are true or false
            %identForm- string specifying whether E is identified or K and
                %G. 'E'=Uses the E formulation. 'KG'= identifies both bulk
                %and shear modulus. NOTE: K and G are always passed into
                %the constitutive model. 
           %minFunc- which function is used to perform the minimization
     %MinOpts- minimization option imputs for the function of choice      
        %Algorithm-  which minimization algorithm is used in the
                %minimization function
     %CondOpts- Data Conditioning options. Spatial Downsampling and data
        %sensoring are perfomed before being passed into the function.
        %Temporal downsampling is performed here. Relevant fields
            %TempDS- indicates whether or not temporal downsampling will be
                %performed with options:
                    %true- Temporal downsampling performed
                    %false- No temporal downsampling
            %Tds- Temporal downsampling factor
     
%output parameters
       %SGparam- parametrs identified with stress gage nonlinear vfm
           %algorithm
       %CostFunc- value of the cost function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Put exact Properties into form used by the cost function
            knownParam=[exactProps.Kinf,exactProps.Ginf,exactProps.nu,....
                SolveOpts.constnu,exactProps.Gi,exactProps.tau];

            %%  Run the minimization algorithm if Einf is assumed to be known
            %% set up lack of constraint
            Aeq=[];
            Beq=[];

            %% run the minimization
            identParam=fmincon(@(constParam)func_ViscoBulkCostSG(SG, ...
                knownParam,constParam,strain,time,CondOpts),...
                guessvec,...
                [],... %A
                [],... %b
                Aeq,... %Aeq
                Beq,... beq
                lbvec,ubvec,... %upper and lower bounds
                [],...
                minOpts);

            CostFunc=func_ViscoBulkCostSG(SG,...
                knownParam,identParam,strain,time,CondOpts);



            SG_param=identParam;

end