function SG_frames=func_GenerateSGMovieFrames(SG,X_vec,time,TestDeg)
%This code is written to generate frames for a video of the propagation of
    %shear and normal stresses calculated with the stresss gauge equations

%Author: Andrew Matejunas
%Date Completed: 2022-06-07

%Function input arguments:
    %SG- Structure containing stresses calculated with the stress gauge
        %equations with dimensions [# or X coordinates, # Times] and fields
            %x- Average normal stresses in the X direction
            %s- Average in-plane shear stresses
    %X_vec- vector of X coordinates in the specimen
    %Time- Stucture of time data
    %TestDeg- test designation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Convert Time vector to microseconds and X coordinates to mm
timeMic=time.vec*10^6;
Xmm=X_vec*10^3;

%% Choose Save Directory
saveDir=uigetdir({},'Choose Save Directory for MovieFrames');

%% Define frame title

frameTitle=strcat(saveDir,'/',TestDeg);

%% Generate plots
figure('units','Centimeters','outerposition',[0 0 27.4 17.2])

for k=1:length(timeMic)
plot(Xmm,squeeze(SG.x(:,k))*10^-6)

xlabel('$$X$$ (mm)','interpreter','latex')
ylabel('$$\sigma{}_{xx}$$ (MPa)','interpreter','latex')
ylim([min(SG.x,[],'all'),max(SG.x,[],'all')]*10^-6)
title(strcat(num2str(timeMic(k)),' $$\mu{}s$$'),'interpreter','latex')

saveas(gcf,strcat(frameTitle,'_SGX_',num2str(timeMic(k)),'us.png'))
end

for k=1:length(timeMic)
plot(Xmm,squeeze(SG.s(:,k))*10^-6)

xlabel('$$X$$ (mm)','interpreter','latex')
ylabel('$$\sigma{}_{xy}$$ (MPa)','interpreter','latex')
title(strcat(num2str(timeMic(k)),' $$\mu{}s$$'),'interpreter','latex')
ylim([min(SG.s,[],'all'),max(SG.s,[],'all')]*10^-6)
saveas(gcf,strcat(frameTitle,'_SGXY_',num2str(timeMic(k)),'us.png'))
end



SG_frames=1;
end

