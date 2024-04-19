function [SG_param,CostFunc]=func_PronyShearSGvfm(strainS,time,ShearSG,...
    exactProps,intGuess,ub,lb,SolveOpts,minOpts,CondOpts)

%This script is written to perform a nonlinear identification on a prony
%series viscoelastic model using the stress gage formulation of the
%nonlinear virtual fields method. This particular script is only intended
%to identify the shear moduli and time constant

%Author: Andrew Matejunas

%Created:2023/01/09
%Validated:

%Change log:
  %2023/01/09- This code was adapted from func_PronySGvfmV5. Simultaneous
    %identification of bulk and shear moduli was removed and now shear
    %modulus and time constant will be identified separately with this
    %code. A separate function will be used in conjunction with this code
    %to identify bulk modulus using the previously identified shear.
 %- Also removed any ability to deal with assumptions of constant
    %poisson's ratio or to identify long-term constitutive parameters
  

% Description
    %This script identifies material parameters by performing a nonlinear
    %minimization of the cost function:
        %phi=[(Stress calculated w/ stress gage) - (stress calculated with
        %the constitutive model)]^2
%Input parameters
    %strainS- shear strain, either measured or output by the FE model
    %ShearSG- stress calculated with the stress gage equation
    %exactProps- structure containing exact material properties. Contents
        %will vary depending on the assumptions made in the formulation of
        %the problem. Common parameters are
            %Kinf- Long term bulk modulus (commonly assumed to be the
                %quasi-static bulk modulus)
            %Ginf-long term shear modulus
            %nu- QS Poisson's ratio
    %intGuess- initial guesses for the parameters that will be identified 
        %with order (order is this way to reshape into a vector)
        %[K1,K2,...,Kn;...
        % G1,G2,...,Gn;...
        % tau1,tau2,...,taun];  
           
     %ub- matrix of upper bounds on the parameters w/ rows [Ki;Gi;Taui]
     %lb- matrix of lower bounds on the parameters w/ columns [Ki;Gi;Taui;]
      
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
            knownParam=[exactProps.Kinf,exactProps.Ginf,...
                exactProps.nu,SolveOpts.constnu,exactProps.Ki];

            
%%  Run the minimization algorithm if Einf is assumed to be known
            %% set up bounding vectors
            ubvec=reshape(ub,[1,numel(ub)]);
            lbvec=reshape(lb,[1,numel(lb)]);
            guessvec=reshape(intGuess,[1,numel(intGuess)]);
            
            %% set up lack of constraint
            Aeq=[];
            Beq=[];
            
            %% run the minimization
     identParam=fmincon(@(constParam)func_ViscoShearCost(ShearSG,...
                knownParam,constParam,strainS,time,CondOpts),...
                guessvec,...
                [],... %A
                [],... %b
                Aeq,... %Aeq
                Beq,... beq
                lbvec,ubvec,... %upper and lower bounds
                [],...
                minOpts);

     CostFunc=func_ViscoShearCost(ShearSG,...
                knownParam,identParam,strainS,time,CondOpts);
    
SG_param=identParam;

end