function phi=func_ViscoKGcostSGV10(xxSG,ShearSG,knownParam,constParam,...
    strain,time,CondOpts,wK,wG)
%This code is intended to generate a weighted cost function used for
%identification of the pairs for a vscoelastic constituitive model

%Author: Andrew Matejunas (andrew.matejunas@gmail.com)
%Date Verified:

%Change Log:
    %V5- 2022/07/14 Added ability to temporally downsample data (note this
        %was changed directly from v2)
    %V8- 2023/03/1: Chaged weighting of cost functions to 75/25
        %sigma11/sigma12
    %V9- version skipped to prevent confusion with minimization function
    %V10- 2023/03/02: Added explicit inputs for the weighting of the
        %sigma_11 and sigma_22 portions of the cost function. If no inputs
        %are provided axial and shear components are equally weighted
        
    
%% Function input arguments
    %xxSG- Average stress at each x coordinate calculated using the normal
        %stress gage equation
    %ShearSG- Average stress at each coordinated calculated with the shear
        %stress gage equation
    %knownParam- A row vector of known parameters in the following order
        %Kinf- long term bulk modulus
        %Ginf- long term shear modulus
        %nu- Poisson's ration 
        %constnu- true if nu is constant  
     %constParam- vector 
        %[K1,G1,tau1,...,Kn,Gn,taun]; if Kinf&Ginf are known,
        %[Kinf,Ginf,K1,G1,tau1,...Kn,Gn,taun]; if Kinf&Ginf are unknown
    %CondOpts- Data Conditioning options. Spatial Downsampling and data
        %sensoring are perfomed before being passed into the function.
        %Temporal downsampling is performed here. Relevant fields
            %TempDS- indicates whether or not temporal downsampling will be
                %performed with options:
                    %true- Temporal downsampling performed
                    %false- No temporal downsampling
            %Tds- Temporal downsampling factor
    %wK- weighting for the axial portion of the cost function (fraction of
        %1)
     %wG- weighting for the pure shear portion of the cost function 
        %(fraction of 1)
%% Function output arguments
    %phi- Weighted cost function

%% Convert Known and Unknown Parameters into inputs for constitutive model

if knownParam(1)~=0
     %% QS moduli are assumed known
       %known Parameters
        MatProps.Kinf=knownParam(1);
        MatProps.Ginf=knownParam(2);

        %unknown Parameters
        constParam=reshape(constParam,[numel(constParam)/3,3]); %a little
        %computationally expensive
    
    
 

        MatProps.Ki=constParam(:,1);
        MatProps.Gi=constParam(:,2);
        MatProps.tau=constParam(:,3);
        MatProps.K0=MatProps.Kinf+sum(MatProps.Ki);
        MatProps.G0=MatProps.Ginf+sum(MatProps.Gi);
    elseif knownParam(1)==0
        MatProps.Kinf=constParam(1);
        MatProps.Ginf=constParam(2);
        
       constParam=constParam(3:end);
       constParam=reshape(constParam,[numel(constParam)/3,3]); %a little
        %computationally expensive
        
        MatProps.Ki=constParam(:,1);
        MatProps.Gi=constParam(:,2);
        MatProps.tau=constParam(:,3);
        MatProps.K0=MatProps.Kinf+sum(MatProps.Ki);
        MatProps.G0=MatProps.Ginf+sum(MatProps.Gi);
end

if knownParam(end)==1
    MatProps.nu=knownParam(3);
else
    MatProps.nu=0;
end

% force a calculation with K and G
MatProps.E0=0;
MatProps.Ei=0;
MatProps.Einf=0;

%% Calculate average stresses for the cost function

StressModel=func_ViscoConstitutiveV6(strain.x,strain.y,strain.s,...
    time,MatProps,0,0,0); %note the zeros are here because I am not working
                          %with exact values here
switch CondOpts.TempDS
    case false %Use every time data point.
        Avxy=StressModel.Avxy;
        Avxx=StressModel.Avxx;
    
    case true %Downsample temporally
        [Avxx,Avxy,xxSG,ShearSG]=func_TempDS(StressModel.Avxx,...
            StressModel.Avxy,xxSG,ShearSG,CondOpts.Tds);
end


%% Calculate cost function for K
phi_K=(xxSG-Avxx).^2;

%% Calculate cost function for G
phi_G=(ShearSG-Avxy).^2;

%% Calculate total cost function

%%weight the total cost function
    %If no inputs for cost function weighting exist set phi_G and phi_G
        %weights equal to one another
if exist('wK','var')==0 
    wK=0.5; %
    if exist('wG','var')==0
        % If wG is not explicitly stated calculate its contribution to 
            %total cost function
        wG=1-wK; 
    end
end
phi=(phi_K*wK+phi_G*wG);
phi=sum(phi,'all');

end