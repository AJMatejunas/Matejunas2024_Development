function hf = func_plotTestSSCurve(paramStr,tStr,X,Y)
% Author: L. Fletcher
% PhotoDyn Research Group
% Date Created: 11/12/2018 
% Date Edited: 11/12/2018

linFit = fit(X,Y,'poly1');
XFit = linspace(min(X),max(X),100);
YFit = linFit.p1*XFit + linFit.p2;
hf = figure;
hold on
plot(X,Y,'-+k')
plot(XFit,YFit,'-b')
hold off
title({tStr,[paramStr,'=',sprintf('%0.2f',linFit.p1/10^9)]})
xlabel('Strain Avg.')
ylabel('Accel. Avg.')
legend('FE','Fit')

end

