function StressModel=sfunc_ViscoConstitutiveV5(strainxx,strainyy,strainxy,...
    time, MatProps,K_guess,G_guess,tau_guess)
%Author: Andrew Matejunas
%Date completed (and verified): 11/16/2020
    

%Current version: 2/2/2020
%Change log
    %2/2/2020- reintegrated parfor loop
    %4/20/2021- V5 attempted to fix the code for differeng k, and g values
        %DOES NOT WORK
    %4/30/2021- sfunc solving in parallel causes errors in large data files
        %on some machines switched code to evaluate constitutive model 
        %in serial

%this function is written to compute the stress in an IBII test on a
% viscoelastic material characterized in a prony series formulation
% (generalized maxwell model)

% This code is written to follow the method used in ABAQUS to calculate
% stress from strain as laid out in the ABAQUS documentation (Analysis 3
% Materials) section 22.7.1 time domain viscoelasticity. 

%function input arguments
    %strainxx- array of normal strains in the x direction 
    %strainnyy- array of normal strains in the y direction
    %strainxy-array of normal strains in the xy direction 
    
    %time- vector containing each time point
    %MatProps- A structure containing the constitutive parameters.
        %including:
            %E0- instantaneous elastic modulus (0 if K&G formulation)
            %Einf- long term elastic modulus   (0 if K&G formulation)
            %nu- Elastic modulus ~=0 if constant ((0 if K&G formulation)
            %K0- Instantaneous bulk modulus
            %Kinf- Long term bulk modulus
            %G0- Instaneous Shear modulus
            %Ginf- Long term shear modulus
            %tau- Vecotor of time constants
            %Ei-Vector of spring constants of maxwell elements for each
                %time constant
            %Ki- vector of bulk moduli for each time constant
            %Gi- vector containing shear modulus for each constant
            
   %Function output arguments
    %StressModel- structure of Stress tensor components calculated through 
                  %the viscoelastic constitutive model contains fields
                  %(xx,xy,yy,hy) hy- is hydrostatic stress. Also exports
                  %hydrostatic strain tensor in the hstrain field
                  %nuR- time dependent poisson's ratio
 
% %% Squeeze strain values from 3D array to 2D array                  
%  if length(time)>1
%     strainxx=squeeze(strainxx);
%     strainyy=squeeze(strainyy);
%     strainxy=squeeze(strainxy);
                      
%% assign guess variables if unassigined
if tau_guess~=0
    clear MatProps.tau
    MatProps.tau=tau_guess;
end

if K_guess~=0
    clear Matprops.Ki
    MatProps.Ki=K_guess;
end

if G_guess~=0
    clear MatPRops.Gi
    MatProps.Gi=G_guess;
end

%% Calculate buluk and shear modulus ratios (if E and nu are specified)

 if MatProps.Einf~=0 && MatProps.E0~=0
 %calculate long term shear and bulk modulus
 Ginf=MatProps.Einf/(2*(MatProps.nu+1));
 Kinf=MatProps.Einf/(3*(1-2*MatProps.nu));
 
 %calculate instantaneous shear and bulk modulus
 G0=MatProps.E0/(2*(MatProps.nu+1));
 K0=MatProps.E0/(2*(MatProps.nu+1));
 
%calculate shear and bulk moduls of each element 
 Gi=MatProps.Ei./(2*(MatProps.nu+1));
 Ki=MatProps.Ei./(2*(MatProps.nu+1));

 %if using a K&G formulation
 else
 G0=MatProps.G0;
 Ginf=MatProps.Ginf;
 Gi=MatProps.Gi;
 
 K0=MatProps.K0;
 Kinf=MatProps.Kinf;
 Ki=MatProps.Ki;
 
 end
 
%% calculate time dependent ratios
gi=Gi/G0;
ki=Ki/K0;


gt=zeros(length(Gi),length(time));
kt=zeros(length(Ki),length(time));

for m=1:length(gi)
    for n=1:length(time)
        gt(m,n)=gi(m)*(1-exp(-time(n)/MatProps.tau(m)));
        kt(m,n)=ki(m)*(1-exp(-time(n)/MatProps.tau(m)));
    end
end


gR=1-sum(gt,1); %note when debugging code, ensure that length(gR)=length(time)
kR=1-sum(kt,1);

%Relaxation shear and bulk moduli as a function of time
GR=gR*G0+Ginf;
KR=kR*K0+Kinf;

%% calculate relaxation Poisson's ratio as a function of time.
%These are needed to calculate the stress tensor components from the
%volumetric stress

nuR=(3*KR-2*GR)./(2*(3*KR+GR));



%this may not be the most efficient way to handle this according to ABAQUS
%numerical implementation

 %% Calculate volumetric strain
    %reshape strain arrays into vectors
 strainxxVec=reshape(strainxx,[numel(strainxx),1]);
 strainyyVec=reshape(strainyy,[numel(strainyy),1]);
 %strainxyVec=reshape(strainxy,[numel(strainxy),1]);
 
 %strain in the Z direction for PL-Stress
 %create a vector for the caclulation of strainzz (repeat nu at each time by
%the number of coordinates
n_coord=length(strainxx(:,1,1))*length(strainxx(1,:,1));

 %If initial inputs are in E and nu and nu is assumed constant
 if MatProps.nu~=0
 strainzzVec=-MatProps.nu/(1-MatProps.nu)*(strainxxVec+strainyyVec);
 %if intial inputs are formulated  in G and K, and/or nu is considered to
 %be rate/time dependent
 else

nuRStretch=repmat(nuR,n_coord,1);
nuRStretch=nuRStretch(:);

%compute normal strain vector in the Z direction.
strainzzVec=-nuRStretch./(1-nuRStretch).*(strainxxVec+strainyyVec);    
 end
 
 
 %preallocate memory
 %strain_volVec=zeros(size(strainxxVec));
 
 %might want to figure out how to do this outside of a loop

 

strain_volVec=strainxxVec+strainyyVec+strainzzVec;
    

 
%% Caclulate the gamma_i and evol_i (epsilon^vol_i from ABAQUS documentation) values

%Reshape volumetric strain into two arrays
    %Vector in the same manner as strainxx and strainyy
strain_vol=reshape(strain_volVec,size(strainxx));
    %array of dimensions [# of spacial coordinates,# of time steps]
strain_volTime=reshape(strain_vol,[n_coord,length(time)]);

%Reshape shear strain into array with dimensions 
    %[# of spacial coordinates,# of time steps]
strainxyTime=reshape(strainxy,[n_coord,length(time)]);
strainxxTime=reshape(strainxx,[n_coord,length(time)]);
strainyyTime=reshape(strainyy,[n_coord,length(time)]);


%preallocate memory for speed
gamma_i=zeros(length(gi),length(strainxyTime(:,1)),length(time));
evol_i=zeros(size(gamma_i));
edevxx_i=zeros(size(gamma_i));
edevyy_i=zeros(size(gamma_i));
tau=MatProps.tau;

%calclate deviatoric strains
strainDevX=strainxxTime-strain_volTime/3;
strainDevY=strainyyTime-strain_volTime/3;


%Calculate creep strains according to abaqus (6.14) documentation
for ig=1:length(gi)
      for it=1:length(time)
%         %need to fix indexing
         gamma_i(ig,:,it)=gi(ig)/tau(ig)*trapz(time(1:it),exp(-time(1:it)/tau(ig)).*strainxyTime(:,it-(1:it)'+1),2);
         edevxx_i(ig,:,it)=gi(ig)/tau(ig)*trapz(time(1:it),exp(-time(1:it)/tau(ig)).*strainDevX(:,it-(1:it)'+1),2);
         edevyy_i(ig,:,it)=gi(ig)/tau(ig)*trapz(time(1:it),exp(-time(1:it)/tau(ig)).*strainDevY(:,it-(1:it)'+1),2);
         evol_i(ig,:,it)=ki(ig)/tau(ig)*trapz(time(1:it),exp(-time(1:it)/tau(ig)).*strain_volTime(:,it-(1:it)'+1),2);

     end
end
%remove array dimensions equal to 1
gamma_sum=squeeze(sum(gamma_i,1));
edevxx_sum=squeeze(sum(edevxx_i,1));
edevyy_sum=squeeze(sum(edevyy_i,1));
evol_sum=squeeze(sum(evol_i,1));

%if only one element fix array dimensions
if length(gamma_sum(1,:))==1
    gamma_sum=gamma_sum';
    evol_sum=evol_sum';
end

    
%% Calculate shear part of the stress tensor (full deviatoric stress tensor)
TAU=G0.*(strainxyTime-gamma_sum);
TAU=reshape(TAU,size(strainxx));

%deviatoric stress tensor components for xx and yy directions
Devxx=2*G0.*(strainDevX-edevxx_sum);
Devxx=reshape(Devxx,size(strainxx));
Devyy=2*G0.*(strainDevY-edevyy_sum);
Devyy=reshape(Devyy,size(strainyy));

%% Calculate the volumetric part of the stress tensor (hydrostatic stress)
P=-K0*(strain_volTime-evol_sum);
P=reshape(P,size(strainxx));

%% Calculate and export stress tensor components
%We already know the shear components because they are calculated as TAU
StressModel.xy=TAU;

%Now we must calculate the normal components
StressModel.xx=Devxx-P;
StressModel.yy=Devyy-P;

StressModel.hy=squeeze(P); %export hydrostatic stress

%% calculate and export Average Stress Components at each x slice
StressModel.Avxx=squeeze(mean(StressModel.xx));
StressModel.Avyy=squeeze(mean(StressModel.yy));
StressModel.Avxy=squeeze(mean(StressModel.xy)); 
StressModel.Avh=squeeze(mean(StressModel.hy));

%StressModel.xx(isnan(StressModel.xx))=0;

%% Export Hydrostatic strain-
StressModel.hstrain=squeeze(strain_vol);

%% Export Poisson's ratio evolution
StressModel.nuR=nuR;
 
end