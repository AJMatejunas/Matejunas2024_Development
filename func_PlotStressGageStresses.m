function SGplots=func_PlotStressGageStresses(SG,strain,X_vec,time,TestDeg)
% This code is written to generate stress-time and stress-strain curves for
    %an image based inertial impact test
    
 %Author: Andrew Matejunas
 %Date complete
    
%Function Input argument
     %SG- Structure containing stresses calculated with the stress gauge
        %equations with dimensions [# or X coordinates, # Times] and fields
            %x- Average normal stresses in the X direction
            %s- Average in-plane shear stresses
     %strain- structure of strains with fields
        %x- normal strains in x direction
        %y- normal strains in y direction
        %s- in-plane shear strains
    %X_vec- vector of X coordinates in the specimen
    %Time- Stucture of time data
    %TestDeg- test designation
    
    
%% chose save directory
saveDir=uigetdir({},...
    'Choose Directory in which to save SG-strain and time plots');

%% plot stress-time curves
timeMic=time.vec*10^6;

%% Divide up X_vec

Xnum=length(X_vec);

%interval of X indexes such that 6 x_locations are investigated
Xsep=floor(Xnum/5);
%% Get indexes
Xindexes=1:Xsep:Xnum;
X_coords=X_vec(Xindexes)*10^3; %coordinates in mm

%% Calculate Average Strains
strain.AvgX=squeeze(mean(strain.x));
strain.AvgS=squeeze(mean(strain.s));

%% Plot normal stress time curves

figure('units','Centimeters','outerposition',[0 0 27.4 17.2])
plot(timeMic,SG.x(Xindexes(1),:)*10^-6)
Legend=cell(length(Xindexes),1);
Legend{1}=strcat(num2str(X_coords(Xindexes(1))),' mm');
hold on
for k=2:length(Xindexes)
    plot(timeMic,SG.x(Xindexes(k),:)*10^-6)
    LegendString=num2str(X_coords(k));
    Legend{k}=strcat(LegendString,' mm');
end

hold off

lgd=legend(Legend);
xlabel('Time ($$\mu{}s$$)','interpreter','latex')
ylabel('$$\sigma{}_{xx}$$ (MPa)','interpreter','latex')

saveas(gcf,strcat(saveDir,'/',TestDeg,'_SGx_time.fig'))
saveas(gcf,strcat(saveDir,'/',TestDeg,'_SGx_time.png'))
saveas(gcf,strcat(saveDir,'/',TestDeg,'_SGx_time.svg'))

%% Plot normal stress strain curves
figure('units','Centimeters','outerposition',[0 0 27.4 17.2])
plot(strain.AvgX(Xindexes(1),:),SG.x(Xindexes(1),:)*10^-6)
Legend=cell(length(Xindexes),1);
Legend{1}=strcat(num2str(X_coords(Xindexes(1))),' mm');
hold on
for k=2:length(Xindexes)
    plot(strain.AvgX(Xindexes(k),:),SG.x(Xindexes(k),:)*10^-6)
    LegendString=num2str(X_coords(k));
    Legend{k}=strcat(LegendString,' mm');
end

hold off

lgd=legend(Legend)
xlabel('$$\varepsilon{}_{xx}$$ ($$\mu{}s$$)','interpreter','latex')
ylabel('$$\sigma{}_{xx}$$ (MPa)','interpreter','latex')

saveas(gcf,strcat(saveDir,'/',TestDeg,'_SGx_strain.fig'))
saveas(gcf,strcat(saveDir,'/',TestDeg,'_SGx_strain.png'))
saveas(gcf,strcat(saveDir,'/',TestDeg,'_SGx_strain.svg'))


%% Plot shear stress time curves

figure('units','Centimeters','outerposition',[0 0 27.4 17.2])
plot(timeMic,SG.s(Xindexes(1),:)*10^-6)
Legend=cell(length(Xindexes),1);
Legend{1}=strcat(num2str(X_coords(Xindexes(1))),' mm');
hold on
for k=2:length(Xindexes)
    plot(timeMic,SG.s(Xindexes(k),:)*10^-6)
    LegendString=num2str(X_coords(k));
    Legend{k}=strcat(LegendString,' mm');
end

hold off

lgd=legend(Legend);
xlabel('Time ($$\mu{}s$$)','interpreter','latex')
ylabel('$$\sigma{}_{xy}$$ (MPa)','interpreter','latex')

saveas(gcf,strcat(saveDir,'/',TestDeg,'_SGxy_time.fig'))
saveas(gcf,strcat(saveDir,'/',TestDeg,'_SGxy_time.png'))
saveas(gcf,strcat(saveDir,'/',TestDeg,'_SGxy_time.svg'))

%% Plot normal stress strain curves
figure('units','Centimeters','outerposition',[0 0 27.4 17.2])
plot(strain.AvgS(Xindexes(1),:),SG.s(Xindexes(1),:)*10^-6)
Legend=cell(length(Xindexes),1);
Legend{1}=strcat(num2str(X_coords(Xindexes(1))),' mm');
hold on
for k=2:length(Xindexes)
    plot(strain.AvgS(Xindexes(k),:),SG.s(Xindexes(k),:)*10^-6)
    LegendString=num2str(X_coords(k));
    Legend{k}=strcat(LegendString,' mm');
end

hold off

lgd=legend(Legend)
xlabel('$$\varepsilon{}_{xy}$$ ($$\mu{}s$$)','interpreter','latex')
ylabel('$$\sigma{}_{xy}$$ (MPa)','interpreter','latex')

saveas(gcf,strcat(saveDir,'/',TestDeg,'_SGxy_strain.fig'))
saveas(gcf,strcat(saveDir,'/',TestDeg,'_SGxy_strain.png'))
saveas(gcf,strcat(saveDir,'/',TestDeg,'_SGxy_strain.svg'))

SGplots=1;
end