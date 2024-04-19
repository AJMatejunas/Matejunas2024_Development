function hf = func_plotInputPulse(time,specimen,accel)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 10/3/2017
%
% Plots the input pulse from the acceleration over the FOV, uses the
% specimen dimensions and material properties to calculate the input
% contact force
% Create a struct of properties to make nicely formatted figures
plotProps = func_initPlotPropsStruct();

force.x = specimen.mass*squeeze(mean(accel.xAvg));

hf = figure;
set(hf,'Position', [plotProps.locX,plotProps.locY,...
        1.25*2*plotProps.sizePerFigX,1.25*1*plotProps.sizePerFigY])
set(hf,'PaperPositionMode','auto')

% Plot the pulse as acceleration vs time
labelStrs.x = 'Time, (\mus)';
labelStrs.y = 'Avg Accel, a_{x} (m/s^2)';
subplot(1,2,1)
hold on
plot(time.vec*10^6,mean(accel.xAvg),'-xb','linewidth',plotProps.lw,'markersize',plotProps.ms)
xlabel(labelStrs.x,'fontsize',plotProps.hfs,'fontname',plotProps.ft)
ylabel(labelStrs.y,'fontsize',plotProps.hfs,'fontname',plotProps.ft)
set(gca, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
set(gca,'XMinorTick','on','YMinorTick','on')
box on
hold off

% Plot the pulse as force vs time
labelStrs.y = 'Force, F_{x}, (N)';
subplot(1,2,2)
hold on
plot(time.vec*10^6,force.x,'-xb','linewidth',plotProps.lw,'markersize',plotProps.ms)
xlabel(labelStrs.x,'fontsize',plotProps.hfs,'fontname',plotProps.ft)
ylabel(labelStrs.y,'fontsize',plotProps.hfs,'fontname',plotProps.ft)
set(gca, 'fontsize', plotProps.fs,'fontname',plotProps.ft)
set(gca,'XMinorTick','on','YMinorTick','on')
box on
hold off


end

