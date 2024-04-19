function [SG_param,phiTot,phi_K,phi_G]=func_PronySGvfmV10(strain,time,SG, ...
    ShearSG,...
    exactProps,intGuess,ub,lb,SolveOpts,minOpts,CondOpts)

%This script is written to perform a nonlinear identification on a prony
%series viscoelastic model using the stress gage formulation of the
%nonlinear virtual fields method

%Author: Andrew Matejunas
%Created:
%Validated:

%Change log:
    %V2- Added constraints to K and G
    %2/8/2020- added capability to attempt to identify quasistatic modulus
    %V3- Added ability to simultaneously identify K and G
    %V4- Changed constitutive stress reconstruction algorithm to V6 which
        %uses the process described in Mun2006 (V3 is obsolete and cannot
        %identify K and G if effective Poisson's ratio is time dependent.
    %V5- 2022-07-14: Added Temporal Downsampling within the cost function
        %2022-10-11: Added output of the cost function value
    %V6- 2022-12-07: Removed ability to assume constant nu and to identify
                %long term constitutive parameters for computational
                %performance
            %Non-dimentionalized constitutive parameters to make
                %all solved parameters the same order of magnitude
                %in an attempt to speed up code (IT was not
                %successful, and reverted)
  %V8- 2023/02/27: Added output of the bulk and shear portions of the cost
                    %function  
  %V9- 2023/03/01: Added ability to weight the cost function for K and G
                    %and changed the cost function algorithm to 
                    %func_ViscoKGcostSGV8
  %V10-2023/03/02: Added inputs for the weighting of the shear and bulk
                    %moudlus components of the cost function to solveOpts
                    %changed cost function version to func_ViscoKGcostSGV10

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
           %wK- weighting of the sigma_11 portion of the cost function
                    %(some fraction of 1)
                    
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
                SolveOpts.constnu];

            
%%  Run the minimization algorithm if Einf is assumed to be known
if SolveOpts.identEinf==false
    if SolveOpts.constnu==1
        if SolveOpts.KGsame==true
            %% Reshape bounds and initial guesses into vectors
            ubvec=reshape(ub,[1,numel(ub)]);
            lbvec=reshape(lb,[1,numel(lb)]);
            guessvec=reshape(intGuess,[1,numel(intGuess)]);
            
           %% Constrain realtionship between K and G
           nu=exactProps.nu;

           Aeq=[-3*(1-2*nu)/(2*(1+nu)),1,0];
           Beq=0;
            %% run minimization
            identParam=fmincon(@(constParam)func_ViscoCostSGV4(SG,...
                knownParam,constParam,strain,time),...
                guessvec,...
                [],... %A
                [],... %b
                Aeq,... %Aeq
                Beq,... beq
                lbvec,ubvec,... %upper and lower bounds
                [],...
                minOpts);
                
        end    
    end
    if SolveOpts.constnu==0 && SolveOpts.KGsame==false
            %% set up bounding vectors
            ubvec=reshape(ub,[1,numel(ub)]);
            lbvec=reshape(lb,[1,numel(lb)]);
            guessvec=reshape(intGuess,[1,numel(intGuess)]);
            
            %% set up lack of constraint
            Aeq=[];
            Beq=[];
            
            %% run the minimization
            if isfield(SolveOpts,'wK')==0
                SolveOpts.wK=0.5;
            end
if isfield(SolveOpts,'wG')==0
    SolveOpts.wG=1-SolveOpts.wK;
end
     identParam=fmincon(@(constParam)func_ViscoKGcostSGV10(SG,ShearSG,...
                knownParam,constParam,strain,time,CondOpts, ...
                SolveOpts.wK,SolveOpts.wG),...
                guessvec,...
                [],... %A
                [],... %b
                Aeq,... %Aeq
                Beq,... beq
                lbvec,ubvec,... %upper and lower bounds
                [],...
                minOpts);

     [phiTot,phi_K,phi_G]=func_EvalViscoKGcostSGV10(SG,ShearSG,...
                knownParam,identParam,strain,time,CondOpts, ...
                SolveOpts.wK,SolveOpts.wG); %weighting for the different 
                    %portions of the cost function

    end
end

%% Ru the minimization algorithm if Einf is unknown
if SolveOpts.identEinf==true
   if SolveOpts.constnu==1
       if SolveOpts.KGsame==true
          
           
           
           %% Reshape bounds and initial guesses into vectors
           ubvec=reshape(ub,[1,numel(ub)]);
           lbvec=reshape(lb,[1,numel(lb)]);
           guessvec=reshape(intGuess,[1,numel(intGuess)]);
           
           %% Constrain realtionship between K and G
           nu=exactProps.nu; 
           Aeq=[-3*(1-2*nu)/(2*(1+nu)),1,-3*(1-2*nu)/(2*(1+nu)),1,0];
           Beq=0;
           
           %%run minimization
           identParam=fmincon(@(constParam)func_ViscoCostSGV3(SG,...
               knownParam,constParam,strain,time),...
               guessvec,...
               [],...
               [],...
               Aeq,...
               Beq,...
               lbvec,ubvec,...
               [],...
               minOpts);
           
       end
   end
end

SG_param=identParam;

end