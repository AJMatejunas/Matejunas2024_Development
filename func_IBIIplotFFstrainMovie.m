function func_IBIIplotFFstrainMovie(pos,time,disp,SmoothStrain,smoothOpts,...
    extrapOpts,ExpDesig,MovieDir)

%This script is written to make full field movies of the strains from an 
    %IBII test 

%Author: Andrew Matejunas

%Version History/Change log:
    %2023-06-24: Original version

%Function input arguments
    %pos- structure containing coordinates
    %time- structure containing time information
    %disp- structure containing displacement components with and without
        %smoothing
    %SmoothStrain- structure containing smoothed and corrected strain
        %fields
    %SmoothOpts- Structure containing smoothing information
    %extrapOpts- Structure containing kinematic field edge extrapolation
        %options
    %ExpDesig- Test designation
    %SavePath- Path to save movie frames

%Function outputs
    %Frames for a movie of Full-field strains
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define number of colormap levels
Nbits=12;
Nlevels=2^Nbits;
%% Create directories
PngDir=strcat(MovieDir,'/strain/PNG');
FigDir=strcat(MovieDir,'/strain/FIG');

mkdir(PngDir);
mkdir(FigDir);

%% Convert pos to mm and time to microsecond
tms=time.vec*10^6; %time vector in microseconds
Xmm=pos.x*10^3; %Vector of X coordinates in mm
Ymm=pos.y*10^3; %Vector of Y coordinates in mm

%% Convert Spatial and Temporal Kernel Strings for Saved figures and labels
SK=smoothOpts.strain.spatialKernelSize(1);
SKs=num2str(SK);

%% Set default settings for figure
    %note that the figure is sized to be played full-screen on a computer
figFont='Arial';
figFsize=18;

labelProps={'Interpreter','latex','FontName',figFont,'FontSize',figFsize};
MajorProps={'Interpreter','latex','FontName',figFont,'FontSize',20};
axprops={'Ydir','Normal','FontSize',18};

FSfig={'Units','Normalized','OuterPosition',[0,0,1,1]};


%% Calculate Strain without smoothing
NoSmoothOpts=smoothOpts;
NoSmoothOpts.strain.spatialSmooth=0;

[NSstrain,~]=func_smoothCalcStrain_v4(pos,time,disp,NoSmoothOpts.strain, ...
    extrapOpts.strain,false);
%% Determine displacment limits
% Sxlim=max(abs(SmoothStrain.x(...
%     extrapOpts.strain.extrapPx1st:(end-extrapOpts.strain.extrapPx1st),...
%     extrapOpts.strain.extrapPx1st:(end-extrapOpts.strain.extrapPx1st),:)),...
%     [],"all")*[-1,1];
Sxlim=[-1,1]*4e-3;
Sylim=max(abs(SmoothStrain.y(...
    extrapOpts.strain.extrapPx1st:(end-extrapOpts.strain.extrapPx1st),...
    extrapOpts.strain.extrapPx1st:(end-extrapOpts.strain.extrapPx1st),:)),...
    [],"all")*[-1,1];
Sslim=max(abs(SmoothStrain.s(...
    extrapOpts.strain.extrapPx1st:(end-extrapOpts.strain.extrapPx1st),...
    extrapOpts.strain.extrapPx1st:(end-extrapOpts.strain.extrapPx1st),:)),...
    [],"all")*[-1,1];


Sxlim=[-1,1]*4e-3;
Sylim=[-1,1]*4e-3;
Sslim=[-1,1]*4e-3;
%% Generate the figure

figure(FSfig{:}) 

for k=1:length(tms)
    Fnum=num2str(k);
    TimeStr=num2str(tms(k));
   
    %Extract strains from the current frame
    NSx=squeeze(NSstrain.x(:,:,k));
    SSx=squeeze(SmoothStrain.x(:,:,k));

    NSy=squeeze(NSstrain.y(:,:,k));
    SSy=squeeze(SmoothStrain.y(:,:,k));

    NSs=squeeze(NSstrain.s(:,:,k));
    SSs=squeeze(SmoothStrain.s(:,:,k));
   
    t=tiledlayout(3,2);
    
    %Strain_xx
    nexttile
    imagesc(Xmm,Ymm,NSx,Sxlim)
    crameri('-roma',Nlevels)
    title('$\varepsilon_{xx}~SK=0,~TK=0$',labelProps{:})
    set(gca,axprops{:})
    colorbar
    

    
    nexttile
    imagesc(Xmm,Ymm,SSx,Sxlim)
    crameri('-roma',Nlevels)
    tstring=strcat('$\varepsilon_{xx}~SK=',SKs,'\mathrm{px}~TK=0$');
    title(tstring,labelProps{:})
    set(gca,axprops{:})
    colorbar
    


    %Strain_yy
    nexttile
    imagesc(Xmm,Ymm,NSy,Sylim)
    crameri('-roma',Nlevels)
    title('$\varepsilon_{yy}~SK=0,~TK=0$',labelProps{:})
    set(gca,axprops{:})
    colorbar
    

    nexttile
    imagesc(Xmm,Ymm,SSy,Sylim)
    crameri('-roma',Nlevels)
    tstring=strcat('$\varepsilon_{yy}~SK=',SKs,'\mathrm{px}~TK=0$');
    title(tstring,labelProps{:})
    set(gca,axprops{:})
    colorbar
    

    %Shear
    nexttile
    imagesc(Xmm,Ymm,NSs,Sslim)
    crameri('-roma',Nlevels)
    title('$\varepsilon_{xy}~SK=0,~TK=0$',labelProps{:})
    title(tstring,labelProps{:})
    set(gca,axprops{:})
    colorbar
    

    nexttile
    imagesc(Xmm,Ymm,SSs,Sslim)
    crameri('-roma',Nlevels)
    tstring=strcat('$\varepsilon_{xy}~SK=',SKs,'\mathrm{px}~TK=0$');
    title(tstring,labelProps{:})
    set(gca,axprops{:})
    colorbar
    

    %Title and axis labels
    Tstring=strcat('Strain Frame~',Fnum,'$~t=',TimeStr,'\mathrm{\mu s}$');
    title(t,Tstring,labelProps{:})
    xlabel(t,'$X \mathrm{(mm)}$',MajorProps{:})
    ylabel(t,'$Y \mathrm{(mm)}$',MajorProps{:})

    PngName=strcat(PngDir,'/',ExpDesig,'_Strain_SK',SKs,'_Frame',Fnum,'.png');
    saveas(gcf,PngName)
    FigName=strcat(FigDir,'/',ExpDesig,'_Strain_SK',SKs,'_Frame',Fnum,'.fig');
    saveas(gcf,FigName)
end

end