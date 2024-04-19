function func_IBBplotFFdisp(pos,time,disp,smoothOpts,ExpDesig,MovieDir)
%This function is written to plot a movie comparing unsmoothed full-field
    %displacements with their smoothed fields from an IBII test

%Author: Andrew Matejunas

%Version History/Change log
    %2023-06-24: Initial version

%Function input arguments
    %pos- structure containing coordinates
    %time- structure containing time information
    %disp- structure containing displacement components with and without
        %smoothing
    %smoothOpts- structure containing smoothing parameters
    %ExpDesig- test desingation
    %MoviDir- Parent directory for all movies

 %Function outputs
    %Frame by frame comparing smoothed and unsmoothed displacements with
        %frames in fig and png format
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create directories
PngDir=strcat(MovieDir,'/disp/PNG');
FigDir=strcat(MovieDir,'/disp/FIG');

mkdir(PngDir);
mkdir(FigDir);

%% Convert pos to mm and tiem to microsecond
tms=time.vec*10^6; %time vector in microseconds
Xmm=pos.x*10^3; %Vector of X coordinates in mm
Ymm=pos.y*10^3; %Vector of Y coordinates in mm

%% Convert Spatial and Temporal Kernel Strings for Saved figures and labels
SK=smoothOpts.strain.spatialKernelSize(1);
SKs=num2str(SK);

TK=smoothOpts.accel.TemporalKernelSize(1);
TKs=num2str(TK);


%% Set default settings for figure
    %note that the figure is sized to be played full-screen on a computer
figFont='Arial';
figFsize=18;

labelProps={'Interpreter','latex','FontName',figFont,'FontSize',figFsize};
MajorProps={'Interpreter','latex','FontName',figFont,'FontSize',20};
axprops={'Ydir','Normal','FontSize',18};

FSfig={'Units','Normalized','OuterPosition',[0,0,1,1]};

%% Determine displacment limits
Dxlim=max(abs(disp.x),[],"all")*[-1,1];
Dylim=max(abs(disp.y),[],"all")*[-1,1];

%% Generate the figure

figure(FSfig)

for k=1:length(tms)
    Fnum=num2str(K);
    TimeStr=num2str(tms(k));
   
    %Extract displacements from the current frame
    Dx=squeeze(disp.x(:,:,k))*10^3;
    DxSS=squeeze(disp.sSmooth.x(:,:,k))*10^3;
    DxTS=squeeze(disp.tSmooth.x(:,:,k))*10^3;

    Dy=squeeze(disp.y(:,:,k))*10^3;
    DySS=squeeze(disp.sSmooth.y(:,:,k))*10^3;
    DyTS=squeeze(disp.tSmooth.y(:,:,k))*10^3;

    t=tiledlayout(3,2);
    
    %No Smoothing
    nexttile
    imagesc(Xmm,Ymm,Dx,Dxlim)
    crameri('vik')
    title('$u_x~\mathrm{(mm)}~SK=0,~TK=0$',labelProps{:})

    nexttile
    imagesc(Xmm,Ymm,Dy,Dylim)
    crameri('vik')
    title('$u_y~\mathrm{(mm)}~SK=0,~TK=0$',labelProps{:})
    
    %Spatial smoothing
    nexttile
    imagesc(Xmm,Ymm,DxSS,Dxlim)
    crameri('vik')
    tstring=strcat('$u_x~\mathrm{(mm)}~SK=',SKs,'\mathrm{px}~TK=0$');
    title(tstring,labelProps{:})

    nexttile
    imagesc(Xmm,Ymm,DySS,Dylim)
    crameri('vik')
    tstring=strcat('$u_y~\mathrm{(mm)}~SK=',SKs,'\mathrm{px}~TK=0$');
    title(tstring,labelProps{:})

    %Temporal Smoothing
    imagesc(Xmm,Ymm,DxTS,Dxlim)
    crameri('vik')
    tstring=strcat('$u_x~\mathrm{(mm)}~SK=0~TK=',TKs,'\mathrm{frames}$');
    title(tstring,labelProps{:})

    nexttile
    imagesc(Xmm,Ymm,DyTS,Dylim)
    crameri('vik')
    tstring=strcat('$u_y~\mathrm{(mm)}~SK=0~TK=',TKs,'\mathrm{frames}$');
    title(tstring,labelProps{:})


    %Title and axis labels
    Tstring=strcat('Displacement Frame~',Fnum,'$~t=',TimeStr,'\mathrm{\mu s}$');
    title(t,Tstring,labelProps{:})
    xlabel(t,'$X \mathrm{(mm)}$',labelProps{:})
    ylabel(t,'$Y \mathrm{(mm)}$',labelProps{:})

    PngName=strcat(PngDir,'/',ExpDesig,'_Displacement_Frame',Fnum,'.png');
    saveas(gcf,PngName)
    FigName=strcat(FigDir,'/',ExpDesig,'_Displacement_Frame',Fnum,'.fig');
    saveas(gcf,FigName)
    
end






end