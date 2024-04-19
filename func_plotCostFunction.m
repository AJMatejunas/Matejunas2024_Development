function [phiKTax,phiKTtot,phiKGax,phiKGtot,phiGTax,phiGTtot,phiGTs] =...
    func_plotCostFunction(SG,knownParam,strain,timevec,CondOpts, ...
    wK,wG,...
    refParams,identParams,ub,lb,SaveDir,Desig,Source)
%This function is written to produce plots of the cost function for a given
    %range of input parameters. 

%Author: Andrew Matejunas
%Date Created: 2023/03/04

%Version History/Change Log:
    %phiKTax- Axial compontent of the cost function varying K and tau

%Function Input arguments:
    %SG- structure containing average stress gauge stresses with fields
        %x- average axial stresses
        %s- average shear stresses
    %knownParam- A row vector of known constitutitve parameters in the
        %following order
            %Kinf- long term bulk modulus
            %Ginf- long term shear modulus
            %nu- Poisson's ration 
            %constnu- true if nu is constant
    %constParam- vector 
        %[K1,G1,tau1,...,Kn,Gn,taun]; if Kinf&Ginf are known,
        %[Kinf,Ginf,K1,G1,tau1,...Kn,Gn,taun]; if Kinf&Ginf are unknown
    %strain- structure containing the strain fields for calculation of the
        %constitutive model stresses (FE output or calculated from the grid
        %method)
            %x- axial strain 
            %y- transverse strain
            %s- shear strain
    %timevec-vector containing each time point 
    %CondOPts- Data Conditioning options. Spatial Downsampling and data
        %censoring are perfomed before being passed into the function.
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
    %refParams- reference constitutive parameters structure with fields
        %K- reference bulk moduli
        %G- reference shear modulus
        %tau- reference time constant
    %identParams- constiutive parameter structure of identfied parameters
        %K- bulk modulus
        %G- shear modulus
        %tau- time constant
     %ub- vector of upper bounds of parameter ident/reference values
        %[K,g,tau]
     %lb- vector of lower bounds of parameter ident/reference values
        %[K,g,tau]
     %Desig- string specifying the test designation for the plots
     %source- string specifying data source with options
        %'FE'- data comes from finite element outputs
        %'GM'- data comes from the grid method

%Function Output arguments:
    %phiKTax- Axial compontent of the cost function with identified G and
        %varying K and tau
    %phiKTtot- weighted total of the cost function with identified G and
        %varying K and tau
    %phiGTax- Axial compontent of the cost function with identified K and
        %varying G and tau
    %phiGTtot- weighted total of the cost function with identified K and
        %varying G and tau
    %phiGTs- Shear compontent of the cost function with identified K and
        %varying G and tau
    %phiKGax- Axial compontent of the cost function with identified tau and
        %varying K and G
    %phiKGtot- weighted total of the cost function with identified tau and
        %varying K and G
    %also outputs plots of all of these cost functions both directly and
        %with the logarithm 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set time to form used in constitutive model calculations
time=timevec;

%% Generate vector of parameter inputs
KRatVec=linspace(lb(1),ub(1),51);
KVec=KRatVec*refParams.K;

GRatVec=linspace(lb(2),ub(2),51);
GVec=GRatVec*refParams.G;

tauRatVec=linspace(lb(3),ub(3),51);
tauVec=tauRatVec*refParams.tau;


%% Pull out Gident for parfor
    G=identParams.G;
    tau=identParams.tau;
    K=identParams.K;
%% Initialize cost functions
%bulk-time constant
phiKTax=zeros(length(KRatVec),length(tauRatVec));
phiKTtot=zeros(size(phiKTax));

%bulk-shear 
phiKGax=zeros(length(KRatVec),length(GRatVec));
phiKGtot=zeros(size(phiKGax));

%shear- time constant
phiGTax=zeros(length(GRatVec),length(tauRatVec));
phiGTtot=zeros(size(phiGTax));
phiGTs=zeros(size(phiGTtot));
%% reomove SG components for use in parfor loop
xxSG=SG.x;
ShearSG=SG.s;
%% Generate cost functions with bulk
for k=1:length(KRatVec)
    tempK=KVec(k);
    parfor t=1:length(tauRatVec)
   constParam=[tempK,G,tauVec(t)];
   [phiKTtot(k,t),phiKTax(k,t),~]=func_EvalViscoKGcostSGV10(xxSG,ShearSG, ...
        knownParam, ...
        constParam,...
        strain,time,CondOpts,wK,wG);
    end
    parfor g=1:length(GRatVec)
    constParam=[tempK,GVec(g),tau];
   [phiKGtot(k,g),phiKGax(k,g),~]=func_EvalViscoKGcostSGV10(xxSG,ShearSG, ...
        knownParam, ...
        constParam,...
        strain,time,CondOpts,wK,wG);
    end 
end

%% Generate G-tau cost functions
for g=1:length(GVec)
    tempG=GVec(g);
    parfor t=1:length(tauRatVec)
        constParam=[K,tempG,tauVec(t)];
        [phiGTtot(g,t),phiGTax(g,t),phiGTs(g,t)]=...
            func_EvalViscoKGcostSGV10(xxSG,ShearSG, ...
        knownParam, ...
        constParam,...
        strain,time,CondOpts,wK,wG);
    end
end

%% calculate logarithms of the cost functions
logPhi.KTtot=log10(phiKTtot);
logPhi.KTax=log10(phiKTax);
logPhi.KGtot=log10(phiKGtot);
logPhi.KGax=log10(phiKGax);
logPhi.GTtot=log10(phiGTtot);
logPhi.GTax=log10(phiGTax);
logPhi.GTs=log10(phiGTs);

%% PLot K-tau cost functions
figure('units','centimeters','InnerPosition',[0,0,18.2,18.2])
    %Total cost fucntion
    subplot(2,2,1)
    contourf(KRatVec,tauRatVec,phiKTtot',50,'LineStyle','none')
    title('\phi(K,\tau)_{total}')
    xlabel('K_1/K_{1,ref}')
    ylabel('\tau_1/\tau_{1,ref}')
    colormap('cool')
    cx=colorbar;
    cx.Label.String='\phi';
    set(gca,'FontSize',14);
    
    %axial only
    subplot(2,2,2)
    contourf(KRatVec,tauRatVec,phiKTax',50,'LineStyle','none')
    title('\phi(K,\tau)_{axial}')
    xlabel('K_1/K_{1,ref}')
    ylabel('\tau_1/\tau_{1,ref}')
    colormap('cool')
    cx=colorbar;
    cx.Label.String='\phi';
    set(gca,'FontSize',14);
    
    %Total log-scale
    subplot(2,2,3)
    contourf(KRatVec,tauRatVec,logPhi.KTtot',50,'LineStyle','none')
    title('\phi(K,\tau)_{total}')
    xlabel('K_1/K_{1,ref}')
    ylabel('\tau_1/\tau_{1,ref}')
    colormap('cool')
    cx=colorbar;
    cx.Label.String='log_{10}(\phi)';
    set(gca,'FontSize',14);
    set(gca,'ColorScale','log')
    
    %axial only log scale
    subplot(2,2,4)
    contourf(KRatVec,tauRatVec,logPhi.KTax',50,'LineStyle','none')
    title('\phi(K,\tau)_{axial}')
    xlabel('K_1/K_{1,ref}')
    ylabel('\tau_1/\tau_{1,ref}')
    colormap('cool')
    cx=colorbar;
    cx.Label.String='log_{10}(\phi)';
    set(gca,'FontSize',14);
    set(gca,'ColorScale','log')

    %Save figure
    figCall='phiKT';
    figSaveName=strcat(SaveDir,'/',Desig,'_',figCall,'_',Source);
    saveas(gcf,figSaveName,'fig')
    saveas(gcf,figSaveName,'png')
    saveas(gcf,figSaveName,'svg')
%% Plot K-G cost functions
figure('units','centimeters','InnerPosition',[0,0,18.2,18.2])
    %Total cost fucntion
    subplot(2,2,1)
    contourf(KRatVec,GRatVec,phiKGtot',50,'LineStyle','none')
    title('\phi(K,G)_{total}')
    xlabel('K_1/K_{1,ref}')
    ylabel('G_1/G_{1,ref}')
    colormap('cool')
    cx=colorbar;
    cx.Label.String='\phi';
    set(gca,'FontSize',14);
    
    %axial only
    subplot(2,2,2)
    contourf(KRatVec,GRatVec,phiKGax',50,'LineStyle','none')
    title('\phi(K,G)_{axial}')
    xlabel('K_1/K_{1,ref}')
    ylabel('G_1/G_{1,ref}')
    colormap('cool')
    cx=colorbar;
    cx.Label.String='\phi';
    set(gca,'FontSize',14);
    
    %Total log-scale
    subplot(2,2,3)
    contourf(KRatVec,GRatVec,logPhi.KGtot',50,'LineStyle','none')
    title('\phi(K,G)_{total}')
    xlabel('K_1/K_{1,ref}')
    ylabel('G_1/G_{1,ref}')
    colormap('cool')
    cx=colorbar;
    cx.Label.String='log_{10}(\phi)';
    set(gca,'FontSize',14);
    set(gca,'ColorScale','log')

    %axial only log scale
    subplot(2,2,4)
    contourf(KRatVec,GRatVec,logPhi.KGax',50,'LineStyle','none')
    title('\phi(K,G)_{axial}')
    xlabel('K_1/K_{1,ref}')
    ylabel('G_1/G_{1,ref}')
    colormap('cool')
    cx=colorbar;
    cx.Label.String='log_{10}(\phi)';
    set(gca,'FontSize',14);
    set(gca,'ColorScale','log')
    
    %Save figure
    figCall='phiKG';
    figSaveName=strcat(SaveDir,'/',Desig,'_',figCall,'_',Source);
    saveas(gcf,figSaveName,'fig')
    saveas(gcf,figSaveName,'png')
    saveas(gcf,figSaveName,'svg')
%% plot G-tau cost functions
figure('units','centimeters','InnerPosition',[0,0,18.2,27.3])
    %Total cost fucntion
    subplot(3,2,1)
    contourf(GRatVec,tauRatVec,phiGTtot',50,'LineStyle','none')
    title('\phi(G,\tau)_{total}')
    xlabel('G_1/G_{1,ref}')
    ylabel('\tau_1/\tau_{1,ref}')
    colormap('cool')
    cx=colorbar;
    cx.Label.String='\phi';
    set(gca,'FontSize',14);
    
    %axial only
    subplot(3,2,3)
    contourf(GRatVec,tauRatVec,phiGTax',50,'LineStyle','none')
    title('\phi(G,\tau)_{axial}')
    xlabel('G_1/G_{1,ref}')
    ylabel('\tau_1/\tau_{1,ref}')
    colormap('cool')
    cx=colorbar;
    cx.Label.String='\phi';
    set(gca,'FontSize',14);

    %Total log-scale
    subplot(3,2,2)
    contourf(GRatVec,tauRatVec,logPhi.GTtot',50,'LineStyle','none')
    title('\phi(G,\tau)_{total}')
    xlabel('G_1/G_{1,ref}')
    ylabel('\tau_1/\tau_{1,ref}')
    colormap('cool')
    cx=colorbar;
    cx.Label.String='log_{10}(\phi)';
    set(gca,'FontSize',14);
    set(gca,'ColorScale','log')
    
    %axial only log scale
    subplot(3,2,4)
    contourf(GRatVec,tauRatVec,logPhi.GTax',50,'LineStyle','none')
    title('\phi(G,\tau)_{axial}')
    xlabel('G_1/G_{1,ref}')
    ylabel('\tau_1/\tau_{1,ref}')
    colormap('cool')
    cx=colorbar;
    cx.Label.String='log_{10}(\phi)';
    set(gca,'FontSize',14);
    set(gca,'ColorScale','log')
    
    %shear only cost fucntion
    subplot(3,2,5)
    contourf(GRatVec,tauRatVec,phiGTs',50,'LineStyle','none')
    title('\phi(G,\tau)_{shear}')
    xlabel('G_1/G_{1,ref}')
    ylabel('\tau_1/\tau_{1,ref}')
    colormap('cool')
    cx=colorbar;
    cx.Label.String='\phi';
    set(gca,'FontSize',14);

    %shear only log
    subplot(3,2,6)
    contourf(GRatVec,tauRatVec,logPhi.GTs',50,'LineStyle','none')
    title('\phi(G,\tau)_{shear}')
    xlabel('G_1/G_{1,ref}')
    ylabel('\tau_1/\tau_{1,ref}')
    colormap('cool')
    cx=colorbar;
    cx.Label.String='log_{10}(\phi)';
    set(gca,'FontSize',14);
    set(gca,'ColorScale','log')

    %Save figure
    figCall='phiGT';
    figSaveName=strcat(SaveDir,'/',Desig,'_',figCall,'_',Source);
    saveas(gcf,figSaveName,'fig')
    saveas(gcf,figSaveName,'png')
    saveas(gcf,figSaveName,'svg')
end