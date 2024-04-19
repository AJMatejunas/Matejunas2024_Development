% This script is written to plot a range of cost functions for a
% viscoelastic material characterized by a generalized Maxwell model with 2
% Maxwell elements

%% Initialize
clear variables; close all; clc

%% Set default fonts for figures to arial for consitency with inkscape
figFont='Arial';
figFsize=10;

labelProps={'Interpreter','latex','FontName',figFont,'FontSize',figFsize};
MajorProps={'Interpreter','latex','FontName',figFont,'FontSize',14};
axprops={'Ydir','Normal','FontSize',10};
%% Set figure size variables
%Set figure widths (Note these come from strain author guidelines
SCW=8.3; %Single Column in cm
DCW=17.5; %Double Column Width in cm
CW15=12.5; %1.5 column width

AR11=1; %standard aspect ratio of 1 to 1
AR43=3/4; %4:3 aspect ration for line plots
AR23=2/3;
AR32=3/2;
%single figure (defaults to 1.5 column width) 
    %ideal for single 1-D plots

SingFig={'Units','Centimeters','InnerPosition',[1,1,CW15,AR43*CW15]};

%Double Colum figure for Double column width
DcSrFig={'Units','Centimeters','InnerPosition',[1,1,DCW,AR11*SCW]};
%Double column double row
DcDrFig={'Units','Centimeters','InnerPosition',[1,1,DCW,AR11*DCW]};
%Triple column
TcDrFig={'Units','Centimeters','InnerPosition',[1,1,DCW,AR23*DCW]};
DcTrFig={'Units','Centimeters','InnerPosition',[1,1,DCW,AR11*DCW]};
T3B2Fig={'Units','Centimeters','InnerPosition',[1,1,DCW,AR32*DCW]};
%Full Page Portrait Figure
FPpFig={'Units','Centimeters','InnerPosition',[1,1,DCW,20]};
%Full Page Landscape
FPlFig={'Units','Centimeters','InnerPosition',[0.5,0.5,230,DCW]};


%% Choose folder in which to save results
SaveDir=uigetdir('','Choose Folder to Save Data');

%% Select files to load
[SGname,SGpath]=uigetfile('*.mat','Choose file containing SG data');
SGfile=strcat(SGpath,'/',SGname);

[PropsName,PropsPath]=uigetfile('*.mat', ...
    'Choose file containing material properties');
Propsfile=strcat(PropsPath,'/',PropsName);

%% load files
fprintf('Loading Stress Gauge Data \n');
load(SGfile);

fprintf('Loading material properties file \n')
load(Propsfile);

%% Put SG and stain structs into form used by the cost function algorithm
SG=Full_SG;
ShearSG=SG.s;
strainS=strain.s;

%% CondOpts Structure
CondOpts.FreeCens=20;
CondOpts.cutStartFrame=0;
CondOpts.cutEndFrames=4;
CondOpts.TempDS=0;

%% Generate Fractional Vectors
FracVec=0.75:0.01:1.25;

K1vec=FracVec*MatProps.Ki(1);
K2vec=FracVec*MatProps.Ki(2);

G1vec=FracVec*MatProps.Gi(1);
G2vec=FracVec*MatProps.Gi(2);

tau1vec=FracVec*MatProps.tau(1);
tau2vec=FracVec*MatProps.tau(2);

%% Convert frational vector to grids for shear modulus

[G2gridy,G1gridx]=meshgrid(G2vec,G1vec);
[tau2gridy,tau1gridx]=meshgrid(tau2vec,tau1vec);

[tau1gridy,G2gridx]=meshgrid(tau1vec,G2vec);


%% Initialize the cost functions
[phiG1G2,phiT1T2,phiG1T1,phiG2T2,phiG1T2,phiG2T1]=...
    deal(zeros(size(G1gridx)));
%% Create structute of known parameters
    Einf=2.98E9;
    nu=0.26;

    exactProps.Kinf=Einf/(3*(1-2*nu));
    exactProps.Ginf=Einf/(2*(1+nu));
    exactProps.nu=0;
    exactProps.Ki=MatProps.Ki;
   
    knownParamShear=[exactProps.Kinf,exactProps.Ginf,...
        exactProps.nu,false,exactProps.Ki];

%% Set array containing ident param

G1G2ParCell=cell(size(G2gridx));
fprintf('Calculating Input Variables for G1-G2 Sweep \n')
parfor k=1:numel(G2gridx)
    G1G2ParCell{k}=[G1gridx(k),MatProps.tau(1);G2gridy(k),MatProps.tau(2)];
end
fprintf('G1-G2 Input Variable calculation complete \n')
%% Calculate the cost function for G1 and G2
fprintf('Calculating G1-G2 CF sweep \n')

parfor k=1:numel(G2gridy)
    constParam=cell2mat(G1G2ParCell(k));
    phiG1G2(k)=func_ViscoShearCost(ShearSG,...
        knownParamShear,constParam,strainS,time.vec,CondOpts);
end

fprintf('G1-G2 Sweep Complete \n')

%% Locate G1G2 minimum cost function
MinG1G2.Value=min(phiG1G2,[],'all');
MinG1G2.ind=find(phiG1G2==MinG1G2.Value);
[MinG1G2.indRow,MinG1G2.indCol]=ind2sub(size(phiG1G2),MinG1G2.ind);

%% Plot G1-G2 Sweep Results
figure(SingFig{:})
imagesc(FracVec,FracVec,phiG1G2)
xlabel('$\frac{G_1}{G_1^\mathrm{ref}}$',labelProps{:});
ylabel('$\frac{G_2}{G_2^\mathrm{ref}}$',labelProps{:});
colorbar
crameri batlowW
hold on
scatter(FracVec(MinG1G2.indCol),FracVec(MinG1G2.indRow),50,'x')
hold off

FigName=strcat(SaveDir,'/G1G2_CostFunction');
saveas(gcf,FigName,'fig')
saveas(gcf,FigName,'svg')
saveas(gcf,FigName,'png')


%% Plot G1-G2 Sweep Results
figure(SingFig{:})
imagesc(FracVec,FracVec,phiG1G2)
xlabel('$\frac{G_1}{G_1^\mathrm{ref}}$',labelProps{:});
ylabel('$\frac{G_2}{G_2^\mathrm{ref}}$',labelProps{:});
cx=colorbar;
crameri batlowW
hold on
scatter(FracVec(MinG1G2.indCol),FracVec(MinG1G2.indRow),50,'x')
hold off
set(gca,'ColorScale','log')

FigName=strcat(SaveDir,'/G1G2_CostFunction_Log');
saveas(gcf,FigName,'fig')
saveas(gcf,FigName,'svg')
saveas(gcf,FigName,'png')

%% Calculate Tau1-Tau2 Prameter structure for tau1-tau2 CF in shear

fprintf('Generating Parameters structure for tau1-tau2 CF in shear \n')
T1T2ParCell=cell(size(tau1gridx));
parfor k=1:numel(T1T2ParCell)
    T1T2ParCell{k}=...
    [MatProps.Gi(1),tau1gridx(k);MatProps.Gi(2),tau2gridy(k)];
end

fprintf('Tau1-Tau2 Shear Params structs complete \n')

%% Calculate the cost function for tau1-tau2 in shear
fprintf('Calculating Cost fucntion for tau1 and tau2 in Shear \n')

parfor k=1:numel(G2gridy)
    constParam=cell2mat(T1T2ParCell(k));
    phiT1T2(k)=func_ViscoShearCost(ShearSG,...
        knownParamShear,constParam,strainS,time.vec,CondOpts);
end

fprintf('Cost Function calculation for Tau1-Tau2 in shear complete')


%% Loacate Tau1-Tau2 Minimum CF in Shear
MinT1T2.Value=min(phiT1T2,[],'all');
MinT1T2.ind=find(phiT1T2==MinT1T2.Value);
[MinT1T2.indRow,MinT1T2.indCol]=ind2sub(size(phiT1T2),MinT1T2.ind);

%% Plot The tau1 Tau2 Cost Function
figure(SingFig{:})
imagesc(FracVec,FracVec,phiG1G2)
xlabel('$frac{\tau_1}{\tau_1^\mathrm{ref}}$',labelProps{:});
ylabel('$frac{\tau_2}{\tau_2^\mathrm{ref}}$',labelProps{:});
colorbar
crameri batlowW
hold on
scatter(FracVec(MinG1G2.indCol),FracVec(MinG1G2.indRow),50,'x')
hold off

FigName=strcat(SaveDir,'/G1G2_CostFunction');
saveas(gcf,FigName,'fig')
saveas(gcf,FigName,'svg')
saveas(gcf,FigName,'png')

%% Log Scale
figure(SingFig{:})
imagesc(FracVec,FracVec,phiG1G2)
xlabel('$frac{\tau_1}{\tau_1^\mathrm{ref}}$',labelProps{:});
ylabel('$frac{\tau_2}{\tau_2^\mathrm{ref}}$',labelProps{:});
cx=colorbar;
crameri batlowW
hold on
scatter(FracVec(MinG1G2.indCol),FracVec(MinG1G2.indRow),50,'x')
hold off
cx.ColorScale='log';

FigName=strcat(SaveDir,'/G1G2_CostFunction_Log');
saveas(gcf,FigName,'fig')
saveas(gcf,FigName,'svg')
saveas(gcf,FigName,'png')