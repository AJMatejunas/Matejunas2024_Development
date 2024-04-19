function [dispTempPad,timeExtrap] = func_temporalPadDisp(...
    extrapOptsA,time,disp,debug)
% Temporally pads or 'extrapolates' displacement data prior to smoothing
% and acceleration calculations.
%
% Author: Lloyd Fletcher
% Image-Based Mechanics Group (IBMG), University of Southampton
% Date Created: 31/7/2020
% Date Edited: 31/7/2020

   % Create some short hand variables that are 'error' free
    pad.pre = round(abs(extrapOptsA.tempPadFrames(1)));
    pad.post = round(abs(extrapOptsA.tempPadFrames(2)));
    pad.tot = pad.pre+pad.post;
    
    [sy,sx,st] = size(disp.x);
    timeExtrap.vec = 0:time.step:(st-1+pad.tot)*time.step;
    timeExtrap.keepRange = 1+pad.pre:length(timeExtrap.vec)-pad.post;
    
    if ~strcmp(extrapOptsA.tempPadMethod,'replicate')    
        pad.preFitF = round(abs(extrapOptsA.tempPadFitWin(1)));
        pad.postFitF = round(abs(extrapOptsA.tempPadFitWin(2))); 

        pst = st+pad.tot;
        dispTempPad.x = nan(sy,sx,st+pad.tot);
        dispTempPad.y = nan(sy,sx,st+pad.tot);
        dispTempPad.x(:,:,1+pad.pre:st+pad.pre) = disp.x;
        dispTempPad.y(:,:,1+pad.pre:st+pad.pre) = disp.y;
        
        fitT.ind1 = 1+pad.pre:pad.pre+pad.preFitF;
        fitT.ind2 = pst-pad.post-pad.postFitF+1:pst-pad.post;
        fitT.tVals1 = timeExtrap.vec(fitT.ind1);
        fitT.tVals2 = timeExtrap.vec(fitT.ind2);
        fitT.extrap1 = timeExtrap.vec(1:pad.pre);
        fitT.extrap2 = timeExtrap.vec(end-pad.post+1:end);
        %figure; imagesc(dispTempPad.x(:,:,fitT.ind2(end))); colorbar; colormap jet; axis image;
        
        if debug
           ppx = round(sx/2);
           ppy = round(sy/2);
           figure; 
           hold on
           plot(timeExtrap.vec*10^6,squeeze(dispTempPad.x(ppy,ppx,:)*10^3))
           plot(timeExtrap.vec*10^6,squeeze(dispTempPad.y(ppy,ppx,:)*10^3))
           hold off
           xlabel('Time [\mus]'); ylabel('Displacement [mm]');
           legend('disp.x','disp.y'); box on; xlim([0,max(timeExtrap.vec)*10^6]);
        end
        
        for xx = 1:sx
            for yy = 1:sy
                fitT.dispX1 = squeeze(dispTempPad.x(yy,xx,fitT.ind1))';
                fitT.dispX2 = squeeze(dispTempPad.x(yy,xx,fitT.ind2))';
                fitT.dispY1 = squeeze(dispTempPad.y(yy,xx,fitT.ind1))';
                fitT.dispY2 = squeeze(dispTempPad.y(yy,xx,fitT.ind2))';
                
                
                if strcmp(extrapOptsA.tempPadMethod,'quadratic')
                    % Pre-pad the displacement fields in time
                    if pad.pre > 0
                        dispTempPad.x(yy,xx,1:pad.pre) = func_quick1DQuadraticExtrap(...
                            fitT.tVals1,fitT.dispX1,fitT.extrap1);
                        dispTempPad.y(yy,xx,1:pad.pre) = func_quick1DQuadraticExtrap(...
                            fitT.tVals1,fitT.dispY1,fitT.extrap1);
                    end
                    
                    % Post-pad the displacement fields in time
                    if pad.post > 0 
                        dispTempPad.x(yy,xx,end-pad.post+1:end) = func_quick1DQuadraticExtrap(...
                            fitT.tVals2,fitT.dispX2,fitT.extrap2);
                        dispTempPad.y(yy,xx,end-pad.post+1:end) = func_quick1DQuadraticExtrap(...
                            fitT.tVals2,fitT.dispY2,fitT.extrap2);
                    end
                else
                    % Pre-pad the displacement fields in time
                    if pad.pre > 0
                        dispTempPad.x(yy,xx,1:pad.pre) = func_quick1DLinearExtrap(...
                            fitT.tVals1,fitT.dispX1,fitT.tVals1(1),fitT.dispX1(1),...
                            fitT.extrap1);
                        dispTempPad.y(yy,xx,1:pad.pre) = func_quick1DLinearExtrap(...
                            fitT.tVals1,fitT.dispY1,fitT.tVals1(1),fitT.dispY1(1),...
                            fitT.extrap1);
                    end
                    
                    % Post-pad the displacement fields in time
                    if pad.post > 0
                        dispTempPad.x(yy,xx,end-pad.post+1:end) = func_quick1DLinearExtrap(...
                            fitT.tVals2,fitT.dispX2,fitT.tVals2(1),fitT.dispX2(1),...
                            fitT.extrap2);
                        dispTempPad.y(yy,xx,end-pad.post+1:end) = func_quick1DLinearExtrap(...
                            fitT.tVals2,fitT.dispY2,fitT.tVals2(1),fitT.dispY2(1),...
                            fitT.extrap2);
                    end

                end
            end
        end
        
    else % Default option is 'replicate'
        dispTempPad.x = disp.x;
        dispTempPad.y = disp.y;
        
        if pad.pre > 0
            dispTempPad.x = padarray(dispTempPad.x,[0,0,pad.pre],...
                extrapOptsA.tempPadMethod,'pre'); 
            dispTempPad.y = padarray(dispTempPad.y,[0,0,pad.pre],...
                extrapOptsA.tempPadMethod,'pre');
        end
        if pad.post > 0
            dispTempPad.x = padarray(dispTempPad.x,[0,0,pad.post],...
                extrapOptsA.tempPadMethod,'post'); 
            dispTempPad.y = padarray(dispTempPad.y,[0,0,pad.post],...
                extrapOptsA.tempPadMethod,'post');
        end
    end
    
    if debug
       ppx = round(sx/2);
       ppy = round(sy/2);
       figure; 
       hold on
       plot(timeExtrap.vec*10^6,squeeze(dispTempPad.x(ppy,ppx,:)*10^3))
       plot(timeExtrap.vec*10^6,squeeze(dispTempPad.y(ppy,ppx,:)*10^3))
       hold off
       xlabel('Time [\mu s]'); ylabel('Displacement [mm]');
       legend('disp.x','disp.y'); box on; xlim([0,max(timeExtrap.vec)*10^6]);
    end 
end

