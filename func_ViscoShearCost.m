function phi=func_ViscoShearCost(ShearSG,knownParam,constParam,...
    StrainXY,time,CondOpts)
% It has been discovered that previous cost functions in which bulk modulus
    % and shear modulus are solved for simultaneously are unacceptably
    % sensitive to noise. This noise sensitivity leads to high errors in
    % identifaction of generalized Maxwell parameter, most notably the bulk
    % modulus. Therefore this script is written to calculate an independent
    % cost function for the shear modulus and time constant.

%Author: Andrew Matejunas (andrew.matejunas@gmail.com)

%Date Verified:

%Change Log:
    %2023/01/09: this code was adapadted from func_ViscoKGcostSGv5
    %2023/03/15: Clarified the input arguments description
    %2023/03/20: Added the ability to discard end frames from the cost
                    %function evaluation. No version change and should
                    %still maintain backwards compatibility 



    
%% Function input arguments
    %ShearSG- Average stress at each coordinated calculated with the shear
        %stress gage equation
    %knownParam- A row vector of known parameters in the following order
        %Kinf- long term bulk modulus
        %Ginf- long term shear modulus
        %nu- Poisson's ration 
        %constnu- true if nu is constant  
     %constParam- matrix
        %[G1,tau1;G2,tau2;....,Gn,taun]; if Kinf&Ginf are known,
     %CondOpts- Data Conditioning options. Spatial Downsampling and data
        %sensoring are perfomed before being passed into the function.
        %Temporal downsampling is performed here. Relevant fields
            %TempDS- indicates whether or not temporal downsampling will be
                %performed with options:
                    %true- Temporal downsampling performed
                    %false- No temporal downsampling
            %Tds- Temporal downsampling factor
            
%% Function output arguments
    %phi- Weighted cost function

%% Convert Known and Unknown Parameters into inputs for constitutive model


     %% QS moduli are assumed known
       %known Parameters
        MatProps.Kinf=knownParam(1);
        MatProps.Ginf=knownParam(2);
        if length(knownParam)==5
            MatProps.Ki=knownParam(5);
        elseif length(knownParam)>5
            MatProps.Ki=knownParam(5:end);
        end
        %unknown Parameters
        constParam=reshape(constParam,[numel(constParam)/2,2]); %a little
        %computationally expensive
        
 

        MatProps.Gi=constParam(:,1);
        MatProps.tau=constParam(:,2);
        MatProps.K0=MatProps.Kinf+sum(MatProps.Ki);
        MatProps.G0=MatProps.Ginf+sum(MatProps.Gi);

if knownParam(4)==1
    MatProps.nu=knownParam(3);
else
    MatProps.nu=0;
end

% force a calculation with K and G
MatProps.E0=0;
MatProps.Ei=0;
MatProps.Einf=0;

%% Calculate average stresses for the cost function

StressModel=func_ViscoShearConstitutive(StrainXY,...
    time,MatProps); %note the zeros are here because I am not working
                          %with exact values here
switch CondOpts.TempDS
    case false %Use every time data point.
        Avxy=StressModel.Avxy;
          
    case true %Downsample temporally
        tempAvxx=ones(size(StressModel.Avxy));
        tempxxSG=ones(size(ShearSG));
        [~,Avxy,~,ShearSG]=func_TempDS(tempAvxx,...
            StressModel.Avxy,tempxxSG,ShearSG,CondOpts.Tds);
end

%% Crop beginning frames from the constitutive model
if isfield(CondOpts,'cutStartFrames')
    if CondOpts.cutStartFrames<=1
    StressModel.Avxy(:,1:CondOpts.cutStartFrames)=[];
    ShearSG(:,1:CondOpts.cutStartFrame)=[];
    end
end
%% Calculate average strains for input into the thresholding algorithm
% strainS=squeeze(mean(StrainXY));
% 
%  %% Threshold data based on the noise floor
%  [ShearSG,Avxy] =...
%     func_thresholdKinFields(ShearSG,Avxy,...
%      strainS,CondOpts);
%% Calculate cost function for K
% phi_K=(xxSG-Avxx).^2;
%% Calculate cost function for G
phi_G=(ShearSG-Avxy).^2;

%% Calculate total cost function
phi=(phi_G);
phi=sum(phi,'all');

end