function func_IBIIplotFFaccelMovie(pos,time,disp,SmoothAccel,smoothOpts,extrapOpts,...
    diffOpts,ExpAesig,MovieAir)
%This script is written to make full field movies of accelerations for an 
    % IBII test. 

%Author: Andrew Matejunas

%Version History/Change log:
    %2023-06-24: Original version

%Function input arguments
    %pos- structure containing coordinates
    %time- structure containing time information
    %accel- structure containing accellacement components with and without
        %smoothing
    %SmoothAccel- Structure containing smoothed and corrected acceleration
        %fields
    %SmoothOpts- Structure containing smoothing information
    %extrapOpts- Structure containing kinematic field edge extrapolation
        %options
    %diffOpts- structure containing temporal differentiation options for
        %calculation of accelerations
    %ExpAesig- Test designation
    %SavePath- Path to save movie frames

%Function outputs
    %Frames for a movie of full field accelerations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Determine number of colormap Levels
Nbits=12;
Nlevels=2^Nbits;
%% Create directories
PngAir=strcat(MovieAir,'/accel/PNG');
FigAir=strcat(MovieAir,'/accel/FIG');

mkdir(PngAir);
mkdir(FigAir);

%% Convert pos to mm and time to microsecond
tms=time.vec*10^6; %time vector in microseconds
Xmm=pos.x*10^3; %Vector of X coordinates in mm
Ymm=pos.y*10^3; %Vector of Y coordinates in mm

%% Convert Spatial and Temporal Kernel Strings for Saved figures and labels
TK=smoothOpts.accel.temporalKernelSize(1);
TKs=num2str(TK);

%% Calculate acceleration without smoothing
NoSmoothOpts=smoothOpts;
NoSmoothOpts.accel.temporalSmooth=false;

[NSaccel,~,~]=func_smoothCalcAccel_v4(pos,time,disp,...
    NoSmoothOpts.accel,extrapOpts.accel,diffOpts,false);
%% Set default settings for figure
    %note that the figure is sized to be played full-screen on a computer
figFont='Arial';
figFsize=18;

labelProps={'Interpreter','latex','FontName',figFont,'FontSize',figFsize};
MajorProps={'Interpreter','latex','FontName',figFont,'FontSize',20};
axprops={'Ydir','Normal','FontSize',18};

FSfig={'Units','Normalized','OuterPosition',[0,0,1,1]};

%% Aetermine accellacment limits
% Axlim=max(abs(SmoothAccel.x),[],"all")*[-1,1];
% Aylim=max(abs(SmoothAccel.y),[],"all")*[-1,1];

Axlim=5e6*[-1,1];
Aylim=5e6*[-1,1];


%% Generate the figure

figure(FSfig{:})

for k=1:length(tms)
    Fnum=num2str(k);
    TimeStr=num2str(tms(k));
   
    %Extract accellacements from the current frame
    Ax=squeeze(NSaccel.x(:,:,k));
    AxTS=squeeze(SmoothAccel.x(:,:,k));

    Ay=squeeze(NSaccel.y(:,:,k));
    AyTS=squeeze(SmoothAccel.y(:,:,k));

    t=tiledlayout(2,2);
    
    %No Smoothing
    nexttile
    imagesc(Xmm,Ymm,Ax,Axlim)
    crameri('-roma',Nlevels)
    title('$a_x~\mathrm{(\frac{m}{s^2})}~SK=0,~TK=0$',labelProps{:})
    set(gca,axprops{:})
    colorbar
    


    nexttile
    imagesc(Xmm,Ymm,Ay,Aylim)
    crameri('-roma',Nlevels)
    title('$a_y~\mathrm{(\frac{m}{s^2})}~SK=0,~TK=0$',labelProps{:})
    set(gca,axprops{:})
    colorbar
    


    %Temporal Smoothing
    nexttile
    imagesc(Xmm,Ymm,AxTS,Axlim)
    crameri('-roma',Nlevels)
    tstring=strcat('$a_x~\mathrm{(\frac{m}{s^2})}~SK=0~TK=',TKs,'\mathrm{frames}$');
    title(tstring,labelProps{:})
    set(gca,axprops{:})
    colorbar

    nexttile
    imagesc(Xmm,Ymm,AyTS,Aylim)
    crameri('-roma',Nlevels)
    tstring=strcat('$a_y~\mathrm{(\frac{m}{s^2})}~SK=0~TK=',TKs,'\mathrm{frames}$');
    title(tstring,labelProps{:})
    set(gca,axprops{:})
    colorbar

    %Title and axis labels
    Tstring=strcat('Acceleration Frame~',Fnum,'$~t=',TimeStr,'\mathrm{\mu s}$');
    title(t,Tstring,labelProps{:})
    xlabel(t,'$X \mathrm{(mm)}$',MajorProps{:})
    ylabel(t,'$Y \mathrm{(mm)}$',MajorProps{:})

    PngName=strcat(PngAir,'/',ExpAesig,'_accel_TK',TKs,'_Frame', ...
        Fnum,'.png');
    saveas(gcf,PngName)
    FigName=strcat(FigAir,'/',ExpAesig,'_accel_TK',TKs, ...
        '_Frame',Fnum,'.fig');
    saveas(gcf,FigName)

end


end