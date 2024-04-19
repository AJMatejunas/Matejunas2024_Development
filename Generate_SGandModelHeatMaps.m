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


%% Choose Save Location
SaveDir=uigetdir({},'Choose Folder to save figures');

%% Load Stress Gage Data
[SGname,SGpath]=uigetfile('*.mat','load file containing stress gage data');
SGfile=strcat(SGpath,'/',SGname);
fprintf('loading stress gage data \n')
SG=load(SGfile);
fprintf('stress gauge data loaded \n')
X_vec=SG.X_vec*1000;
Y_vec=SG.Y_vec*1000;

%% Load FE stress Data
[FEname,FEpath]=uigetfile('*.mat',...
    'Select File Containing Finite Element Data');
FEfile=strcat(FEpath,'/',FEname);
fprintf('loading finite element stress data \n')
FE=load(FEfile,'stress','time','strain');
fprintf('Finite Element Stresses Loaded \n')

%% Load material properties data file
[refName,refPath]=uigetfile('*.mat', ...
    'Choose file containing the reference constitutive parameters');
load(strcat(refPath,'/',refName),'MatProps');

%% Calculate Constitutive Model Data
fprintf('Calculating constitutive model data \n')
Const.StressModel=func_ViscoConstitutiveV6(FE.strain.x,FE.strain.y, ...
    FE.strain.s,FE.time.vec,MatProps,0,0,0);

%% Calculate Error X-t plots
fprintf('calculating average FE stresses \n')
FE.stress.xAv=squeeze(mean(FE.stress.x));
FE.stress.sAv=squeeze(mean(FE.stress.s));

fprintf('calculating SG-FE errors')
SG.Error.x=(FE.stress.xAv-SG.Full_SG.x)./FE.stress.xAv*100;
SG.Error.x(isnan(SG.Error.x)==1)=0;
SG.Error.s=(FE.stress.sAv-SG.Full_SG.s)./FE.stress.sAv*100;
SG.Error.s(isnan(SG.Error.s)==1)=0;


%% Calculate SG_FE normalized error
SG.NormErr.x=(FE.stress.xAv-SG.Full_SG.x)./max(FE.stress.xAv,[],'all')*100;
SG.NormErr.s=(FE.stress.sAv-SG.Full_SG.s)./max(FE.stress.sAv,[],'all')*100;

%% Calculate Error Statistics
SG.NormErr.xMed=median(SG.NormErr.x,'all');
SG.NormErr.xMean=mean(SG.NormErr.x,'all');
SG.NormErr.xStdev=std(SG.NormErr.x,0,'all');


SG.NormErr.sMed=median(SG.NormErr.s,'all');
SG.NormErr.sMean=mean(SG.NormErr.s,'all');
SG.NormErr.xStdev=std(SG.NormErr.s,0,'all');

%% Calculate FE-Constitutive Model Error

FE.stress.yAv=squeeze(mean(FE.stress.y));
Const.Error.xAv=(Const.StressModel.Avxx-FE.stress.xAv)./FE.stress.xAv*100;
Const.Error.YAv=(Const.StressModel.Avyy-FE.stress.yAv)./FE.stress.yAv*100;
Const.Error.sAv=(Const.StressModel.Avxy-FE.stress.sAv)./FE.stress.sAv*100;

Const.MnormErr.xAv=(Const.StressModel.Avxx-FE.stress.xAv)...
    ./max(FE.stress.xAv,[],'all')*100;
Const.MnormErr.yAv=(Const.StressModel.Avyy-FE.stress.yAv)...
    ./max(FE.stress.yAv,[],'all')*100;
Const.MnormErr.sAv=(Const.StressModel.Avxy-FE.stress.sAv)...
    ./max(FE.stress.sAv,[],'all')*100;

%% Plot figure with superimposed stress magnitude, SG Error, and
    %Constitutive model error
timemic=FE.time.vec*10^6;
% Stress Magnitude
figure(FPpFig{:})
subplot(3,2,1)
imagesc(X_vec,timemic,FE.stress.xAv'*10^-6, ...
    max(abs(FE.stress.xAv*10^-6),[],'all')*[-1,1])
ylabel('time $(\mathrm{\mu{}s)}$',labelProps{:})
title('(a)~$\overline{\sigma{}_{11}^{\mathrm{FE}}}^{x_2}$ (MPa)',labelProps{:})
%caxis(max(abs(FE.stress.xAv*10^-6),[],'all')*[-1,1]);
colorbar
crameri vik
set(gca,'Ydir','normal','FontName','Arial','FontSize',10)

subplot(3,2,2)
imagesc(X_vec,timemic,FE.stress.sAv'*10^-6, ...
    max(abs(FE.stress.sAv*10^-6),[],'all')*[-1,1])
%ylabel('time $(\mathrm{\mu{}s)}$',labelProps{:})
title('(b)~$\overline{\sigma{}_{12}^{\mathrm{FE}}}^{x_2}$ (MPa)',labelProps{:})
%caxis(max(abs(FE.stress.xAv*10^-6),[],'all')*[-1,1]);
colorbar
crameri vik
set(gca,'Ydir','normal','FontName','Arial','FontSize',10)

% Stress Gauge Error
subplot(3,2,3)
imagesc(X_vec,timemic,SG.NormErr.x',1.75*[-1,1])
ylabel('time $(\mathrm{\mu{}s)}$',labelProps{:})
title('(c)~$||\mathrm{Error}(\sigma_{11}^\mathrm{SG})||_{\mathrm{max}}$ (\%)', ...
    labelProps{:})
colorbar
crameri vik
set(gca,'Ydir','normal','FontName','Arial','FontSize',10)

subplot(3,2,4)
imagesc(X_vec,timemic,SG.NormErr.s',1.5*[-1,1])
%ylabel('time $(\mathrm{\mu{}s)}$',labelProps{:})
title('(d)~$||\mathrm{Error}(\sigma_{12}^\mathrm{SG})||_{\mathrm{max}}$ (\%)', ...
    labelProps{:})
colorbar
crameri vik
set(gca,'Ydir','normal','FontName','Arial','FontSize',10)

%Constitutive model error
subplot(3,2,5)
imagesc(X_vec,timemic,Const.MnormErr.xAv',0.1*[-1,1])
ylabel('time $(\mathrm{\mu{}s)}$',labelProps{:})
xlabel('$x_0$~(mm)',labelProps{:})
title('(e)~$||\mathrm{Error}(\sigma_{11}^\mathrm{Model})||_{\mathrm{max}}$ (\%)', ...
    labelProps{:})
colorbar
crameri vik
set(gca,'Ydir','normal','FontName','Arial','FontSize',10)

subplot(3,2,6)
imagesc(X_vec,timemic,Const.MnormErr.sAv',0.025*[-1,1])
%ylabel('time $(\mathrm{\mu{}s)}$',labelProps{:})
title('(f)~$||\mathrm{Error}(\sigma_{12}^\mathrm{Model})||_{\mathrm{max}}$ (\%)', ...
    labelProps{:})
xlabel('$x_0$~(mm)',labelProps{:})
colorbar
crameri vik
set(gca,'Ydir','normal','FontName','Arial','FontSize',10)

figName=strcat(SaveDir,'/X_tMagErr_V4');
saveas(gcf,figName,'fig')
saveas(gcf,figName,'eps')
saveas(gcf,figName,'svg')
