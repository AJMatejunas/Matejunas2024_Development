function [Phi,Phi_K,Phi_G] = func_PlotPhiConstTau(K_vec,G_vec,tau,...
    Full_SG,knownParam,strain,time,CondOpts,ParentDesig,SaveDir, ...
    K_exact,G_exact)
%This script is written to calculate and plot the cost function components
    %for a viscoelastic material with a constant time constant input
        %NOTE: THIS CODE ONLY WORKS FOR GENERALIZED MAXWELL MODELS WITH
        %KNOWN LONG TERM PROPERTIES AND A SINGLE MAXWELL ELEMENT!!!!!!!!

    %Author: Andrew Matejunas


    %Date created: 2022/11/03

    %Chagelog/Version History

    %Function input arguments
        %K_vec- Vector of bulk moduli
        %G_vec- Vector os shear moduli
        %Tau- Time constant input
        %Full SG- structure containing the stress gauge stresses with
            %fields
                %x- normal stresses in the x-direction size [#X points, #
                    %time steps]
                %s- avearaged shear stresses [#X points, # time steps]
        %knownParam- Vector of constitutive parameters that are assumed
            %known
                %[Kinf,Ginf,nu,constnu (1 if constant,0 if assumed rate
                %dependent)]
        %strain- structure containing strain components
        %time- structure containing time data
        %CondOpts- data conditioning options
        %ParentDesig- Test Designation
        %SaveDir- directory where figures will be saved
        %K_exact- reference parameter for bulk modulus
        %G_exact- reference parameter for shear modulus


    %Function output arguments
        %Phi- total cost function
        %Phi_K- bulk modulus portion of the cost function (uses normal
            %stresses in the X direction)
        %Phi_G- Shear modulus protion of the cost function (only depends on
            %shear stresses)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Convert time structure to vector form
timevec=time.vec;

%% Convert SG structure to arrays for cost function algorithm
xxSG=Full_SG.x;
ShearSG=Full_SG.s;
%% Transpose Shear modulus vector for easier parfor loop integration
G_vec=G_vec';
%% Create Parameter vectors for Parfor loop integration
Kin=K_vec.*ones(size(G_vec));
Kvec_Full=reshape(Kin,[numel(Kin),1]);
Gin=G_vec.*ones(size(K_vec));
Gvec_Full=reshape(Gin,[numel(Gin),1]);

%% Pre-allocate memory for cost function
Phi_vec=zeros(size(Kvec_Full));
Phi_Kvec=zeros(size(Kvec_Full));
Phi_Gvec=zeros(size(Kvec_Full));

%% Create Grid inputs for the cost function
%Use this line for testing and debugging use parfor for c
%for n=1:numel(Kvec_Full)
parfor n=1:numel(Kvec_Full)
constParam=[Kvec_Full(n),Gvec_Full(n),tau];
[tPhi,tPhi_K,tPhi_G]=func_ViscoKGcostSGV6(xxSG,ShearSG,knownParam, ...
    constParam,strain,timevec,CondOpts);
Phi_vec(n)=tPhi;
Phi_Kvec(n)=tPhi_K;
Phi_Gvec(n)=tPhi_G;
end

Phi=reshape(Phi_vec,size(Kin));
Phi_K=reshape(Phi_Kvec,size(Kin));
Phi_G=reshape(Phi_Gvec,size(Kin));
%% Plot total cost function
figure('units','centimeters','InnerPosition',[10,10,18,9])

colorscheme='jet';
%Linear Scale
subplot(1,2,1)
K_rat=K_vec/K_exact;
G_rat=G_vec/G_exact;
Phi_log=log10(Phi);
contourf(K_rat,G_rat,Phi,50)
xlabel('K_1/K_{1,ref}')
ylabel('G_1,G_{1,ref}')
colormap(colorscheme)
Cscale=colorbar;
Cscale.Label.String='\phi';

%logarithmic Scale
subplot(1,2,2)
contourf(K_rat,G_rat,Phi_log,50)
xlabel('K_1/K_{1,ref}')
ylabel('G_1,G_{1,ref}')
Cscale=colorbar;
Cscale.Label.String='\phi';
set(gca,'ColorScale','log')

%Save
TotSaveName=strcat(SaveDir,'/',ParentDesig,'_KG_Phi');
saveas(gcf,strcat(TotSaveName,'.fig'))
saveas(gcf,strcat(TotSaveName,'.svg'))
saveas(gcf,strcat(TotSaveName,'.png'))

%% Plot K component
figure('units','centimeters','InnerPosition',[10,10,18,9])

subplot(1,2,1)

%Linear Scale
PhiK_log=log10(Phi_K);
contourf(K_rat,G_rat,Phi_K,50)
xlabel('K_1/K_{1,ref}')
ylabel('G_1,G_{1,ref}')
colormap(colorscheme)
Cscale=colorbar;
Cscale.Label.String='\phi_K';

%logarithmic Scale
subplot(1,2,2)
contourf(K_rat,G_rat,PhiK_log,50)
xlabel('K_1/K_{1,ref}')
ylabel('G_1,G_{1,ref}')
Cscale=colorbar;
Cscale.Label.String='\phi_K';
set(gca,'ColorScale','log')

%Save
KSaveName=strcat(SaveDir,'/',ParentDesig,'_KG_BulkPhi');
saveas(gcf,strcat(KSaveName,'.fig'))
saveas(gcf,strcat(KSaveName,'.svg'))
saveas(gcf,strcat(KSaveName,'.png'))

%% Plot G component
figure('units','centimeters','InnerPosition',[10,10,18,9])

subplot(1,2,1)

%Linear Scale
PhiG_log=log10(Phi_G);
contourf(K_rat,G_rat,Phi_G,50)
xlabel('K_1/K_{1,ref}')
ylabel('G_1,G_{1,ref}')
colormap(colorscheme)
Cscale=colorbar;
Cscale.Label.String='\phi_K';

%logarithmic Scale
subplot(1,2,2)
contourf(K_rat,G_rat,PhiG_log,50)
xlabel('K_1/K_{1,ref}')
ylabel('G_1,G_{1,ref}')
Cscale=colorbar;
Cscale.Label.String='\phi_G';
set(gca,'ColorScale','log')

%Save
GSaveName=strcat(SaveDir,'/',ParentDesig,'_KG_ShearPhi');
saveas(gcf,strcat(GSaveName,'.fig'))
saveas(gcf,strcat(GSaveName,'.svg'))
saveas(gcf,strcat(GSaveName,'.png'))

end