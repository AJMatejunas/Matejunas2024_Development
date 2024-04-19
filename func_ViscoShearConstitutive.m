function StressModel=func_ViscoShearConstitutive(strainxy,...
    time, MatProps)
%Author: Andrew Matejunas
%Date completed: 2023/01/09
    

%Current version: 

%Change log
    %2023/01/09- this code is written from the 2022/12/05 version of
    %func_ViscoConstitutiveV6
    %2023/01/11- fixed an issue resulting in extraneous loop iterations and
        %unacceptable computational time.

%this function is written to independently compute the shear stress in an 
% IBII test on a viscoelastic material characterized in a prony series 
% formulation (generalized maxwell model)

%function input arguments
    %strainxy-array of shear strains in the xy direction 
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
    %StressModel- structure containing the calculated in-plane shear
    %stresses
 
% %% Squeeze strain values from 3D array to 2D array                  
%  if length(time)>1
%     strainxx=squeeze(strainxx);
%     strainyy=squeeze(strainyy);
%     strainxy=squeeze(strainxy);

% number of spatial coordinates
n_coord=length(strainxy(:,1,1))*length(strainxy(1,:,1));
                      
% %% assign guess variables if unassigined
% if tau_guess~=0
%     clear MatProps.tau
%     MatProps.tau=tau_guess;
% end
% 
% if K_guess~=0
%     clear Matprops.Ki
%     MatProps.Ki=K_guess;
% end
% 
% if G_guess~=0
%     clear MatPRops.Gi
%     MatProps.Gi=G_guess;
% end

%% Calculate buluk and shear modulus ratios (if E and nu are specified)

 K0=MatProps.K0;
 Ki=MatProps.Ki;
 G0=MatProps.G0;
 Gi=MatProps.Gi;

tau=MatProps.tau;
Num_mx=length(tau); %number of Maxwell elements
%% Calculate Instanteous Stiffness Tensor
C_Visco0=func_ElasKG(K0,G0);

%% Calculate Vicous Stiffness tensors
C_Viscom=func_ElasKG(Ki,Gi);

% %% Calculate out of plane strain
% %fprintf('Calculating strain33 \n')
% [strain33,eps11m,eps22m,eps33m]=func_ViscoIBIIstrain33V2(strainxx,strainyy,... input strains
%     time,C_Visco0,C_Viscom,tau);
    

%  %% Calculate volumetric strain
% %fprintf('Calculating Volumetric Strain \n')
% strain_vol=strain33+strainxx+strainyy;
% 
% %strain_volTime=reshape(strain_vol,[n_coord,length(time)]);

%% Reshape strain arrays for more efficient calculations
    %[# of spacial coordinates,# of time steps]
strain12Time=reshape(strainxy,[n_coord,length(time)]);
% strain11Time=reshape(strainxx,[n_coord,length(time)]);
% strain22Time=reshape(strainyy,[n_coord,length(time)]);
% strain33Time=reshape(strain33,[n_coord,length(time)]);
%% initialize for speed
Stress12Time=zeros(size(strain12Time));

%% Define necessary variables
Delta_t=zeros(length(time),1);
 
 for k=2:length(time)
     Delta_t(k)=time(k)-time(k-1);
 end
%   [B11m,R11m,B22m,B33m,...
%       R22m,R33m,Q11m,Q22m,Q33m,...
%       R12m,B12m]=deal(zeros(n_coord,Num_mx));
[R12m,B12m]=deal(zeros(n_coord,Num_mx));
%% Calculate normal stresses 
%initialize eps_klm
% [eps11m,eps22m,eps33m,eps12m,eps13m,eps23m]...
%     =deal(zeros([n_coord,Num_mx,length(time)]));
[eps12m,eps13m,eps23m]...
    =deal(zeros([n_coord,Num_mx,length(time)]));
% for n=2:(length(time)-1)
% % following commented out for computational efficiency purposes
% 
% 
%     %define strain inputs for epsklm function
% %    strain11input=strain11Time(:,(n-1):n);
% %    strain22input=strain22Time(:,(n-1):n);
% %    strain33input=strain33Time(:,(n-1):n);
%    
% 
% %    [eps11m(:,:,n),eps22m(:,:,n),eps33m(:,:,n)]...
% %        =func_Muneps_mV2(strain11input,strain22input,...
% %        strain33input,... strain inputs for t_{n-1} and t_{n}
% %    eps11m(:,:,n-1),eps22m(:,:,n-1),eps33m(:,:,n-1),... eps_klm(t_{n-1})
% %    Delta_t(n),tau);
%   
% % for m=1:length(tau)
% %         
% %    B11m(:,m)=strain11Time(:,n)...
% %         +exp(-Delta_t(n+1)/tau(m)).*(eps11m(:,m,n)-strain11Time(:,n))...
% %         +(strain11Time(:,n+1)-strain11Time(:,n)).*(1-tau(m)/Delta_t(n+1)...
% %         .*(1-exp(-Delta_t(n+1)/tau(m))));
% %     R11m(:,m)=C_Viscom(1,1,1,1,m)*B11m(:,m); %for stress11 calculation
% %     Q11m(:,m)=C_Viscom(2,2,1,1,m)*B11m(:,m); %for stress22 calculation
% %     
% %     %22 direction
% %    B22m(:,m)=strain22Time(:,n)...
% %         +exp(-Delta_t(n+1)/tau(m)).*(eps22m(:,m,n)-strain22Time(:,n))...
% %         +(strain22Time(:,n+1)-strain22Time(:,n)).*(1-tau(m)/Delta_t(n+1)...
% %         .*(1-exp(-Delta_t(n+1)/tau(m))));
% %     R22m(:,m)=C_Viscom(1,1,2,2,m)*B22m(:,m);
% %     Q22m(:,m)=C_Viscom(2,2,2,2,m)*B22m(:,m);
% %     
% %     %33 direction
% %     B33m(:,m)=strain33Time(:,n)...
% %         +exp(-Delta_t(n+1)/tau(m)).*(eps33m(:,m,n)-strain33Time(:,n))...
% %         +(strain33Time(:,n+1)-strain33Time(:,n)).*(1-tau(m)/Delta_t(n+1)...
% %         .*(1-exp(-Delta_t(n+1)/tau(m))));
% %     R33m(:,m)=C_Viscom(1,1,3,3,m)*B33m(:,m);
% %     Q33m(:,m)=C_Viscom(2,2,3,3,m)*B33m(:,m);
% %     
% %     
% %  end
% % R11=sum(R11m,2);
% % R22=sum(R22m,2);
% % R33=sum(R33m,2);
% % R=R11+R22+R33;
% % 
% % 
% % Q11=sum(Q11m,2);
% % Q22=sum(Q22m,2);
% % Q33=sum(Q33m,2);
% % 
% % Q=Q11+Q22+Q33;
% % 
% % Stress11Time(:,n+1)=C_Visco0(1,1,1,1)*strain11Time(:,n+1)...
% %     +C_Visco0(1,1,2,2)*strain22Time(:,n+1)...
% %     +C_Visco0(1,1,3,3)*strain33Time(:,n+1)-R;
% % 
% % Stress22Time(:,n+1)=C_Visco0(2,2,1,1)*strain11Time(:,n+1)...
% %     +C_Visco0(2,2,2,2)*strain22Time(:,n+1)...
% %     +C_Visco0(2,2,3,3)*strain33Time(:,n+1)-Q;
%  end
% 

%% Calculate shear stresses
for n=2:(length(time)-1)
    %define strain inputs for epsklm function
   strain12input=strain12Time(:,(n-1):n);
     
   [eps12m(:,:,n),eps13m(:,:,n),eps23m(:,:,n)]...
       =func_Muneps_mshear(strain12input,...
   eps12m(:,:,n-1),eps23m(:,:,n-1),eps23m(:,:,n-1),... eps_klm(t_{n-1})
   Delta_t(n),tau);

    for m=1:length(tau)
        
   B12m(:,m)=strain12Time(:,n)...
        +exp(-Delta_t(n+1)/tau(m)).*(eps12m(:,m,n)-strain12Time(:,n))...
        +(strain12Time(:,n+1)-strain12Time(:,n)).*(1-tau(m)/Delta_t(n+1)...
        .*(1-exp(-Delta_t(n+1)/tau(m))));
   R12m(:,m)=C_Viscom(1,2,1,2,m)*B12m(:,m); %for stress11 calculation
     
    end
R12=sum(R12m,2);

Stress12Time(:,n+1)=C_Visco0(1,2,1,2)*strain12Time(:,n+1)...
    -R12;

end


%% Reshape into [number Ycoord, number X coord, time]
StressModel.xy=reshape(Stress12Time,size(strainxy));

%% calculate and export Average Stress Components at each x slice
StressModel.Avxy=squeeze(mean(StressModel.xy)); 

end