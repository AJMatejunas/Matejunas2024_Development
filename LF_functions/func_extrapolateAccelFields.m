function accelExtrap = func_extrapolateAccelFields(...
    subExtrapOpts,pos,accelCrop,posCropX,posCropY,debug)
% Extrapolates the missing edge data in the displacement fields for the
% grid method or DIC.
%
% Author: Lloyd Fletcher
% Image-Based Mechanics Group (IBMG), University of Southampton
% Date Created: 31/7/2020
% Date Edited: 31/7/2020
    
    % Turn off debug mode by default.
    if nargin < 6
        debug = false;
    end
    
    % Get the sizes of the fields for later calculations
    [ax_sy,ax_sx,st] = size(accelCrop.x);
    [ay_sy,ay_sx,~] = size(accelCrop.y);
    
    % Define extrapolation co-ordinates - X displacement field
    posExtrapX.x = min(posCropX.x)-subExtrapOpts.extrapPx1st*pos.xStep:pos.xStep:...
        max(posCropX.x)+subExtrapOpts.extrapPx1st*pos.xStep;
    posExtrapX.y = min(posCropX.y)-subExtrapOpts.extrapPx2nd*pos.yStep:pos.yStep:...
        max(posCropX.y)+subExtrapOpts.extrapPx2nd*pos.yStep;
    % Define extrapolation co-ordinates - Y displacement field
    posExtrapY.x = min(posCropY.x)-subExtrapOpts.extrapPx2nd*pos.xStep:pos.xStep:...
        max(posCropY.x)+subExtrapOpts.extrapPx2nd*pos.xStep;
    posExtrapY.y = min(posCropY.y)-subExtrapOpts.extrapPx1st*pos.yStep:pos.yStep:...
        max(posCropY.y)+subExtrapOpts.extrapPx1st*pos.yStep;
    
    % Initiliase 'strainExtrap' vars for speed
    if strcmp(subExtrapOpts.fieldOpt,'XY')
        % Here we only extrapolate X in the X direction and same for Y, S
        % we extrapolate both!
        accelExtrap.x = nan(ax_sy,ax_sx+2*subExtrapOpts.extrapPx1st,st);
        accelExtrap.x(:,subExtrapOpts.extrapPx1st+1:end-subExtrapOpts.extrapPx1st,:) = accelCrop.x;
        accelExtrap.y = nan(ay_sy+2*subExtrapOpts.extrapPx1st,ay_sx,st);
        accelExtrap.y(subExtrapOpts.extrapPx1st+1:end-subExtrapOpts.extrapPx1st,:,:) = accelCrop.y;
    else
        % Here we extrapolate in both directions, creates nans on all borders
        accelExtrap.x = nan(ax_sy+2*subExtrapOpts.extrapPx2nd,...
            ax_sx+2*subExtrapOpts.extrapPx1st,st);
        accelExtrap.x(1+subExtrapOpts.extrapPx2nd:end-subExtrapOpts.extrapPx2nd,...
            1+subExtrapOpts.extrapPx1st:end-subExtrapOpts.extrapPx1st,:) = accelCrop.x;

        accelExtrap.y = nan(ay_sy+2*subExtrapOpts.extrapPx1st,...
            ay_sx+2*subExtrapOpts.extrapPx2nd,st);
        accelExtrap.y(1+subExtrapOpts.extrapPx1st:end-subExtrapOpts.extrapPx1st,...
            1+subExtrapOpts.extrapPx2nd:end-subExtrapOpts.extrapPx2nd,:) = accelCrop.y;
    end
    [ae_sy,ae_sx,ae_st] = size(accelExtrap.x);
    
    if debug
    dbf = round(ae_st/2);    
    figure; imagesc(flipud(accelExtrap.x(:,:,dbf))); title('accelExtrap.x');
    axis image; colorbar; colormap jet; set(gca,'YDir','normal');
    figure; imagesc(flipud(accelExtrap.y(:,:,dbf))); title('accelExtrap.y');
    axis image; colorbar; colormap jet; set(gca,'YDir','normal');
    end
    
    %--------------------------------------------------------------------------
    % X DIRECTION
    % Fit data structure for fitting X fields in the X direction
    xFitX.Ind1 = 1+subExtrapOpts.extrapPx1st:...
        subExtrapOpts.extrapPx1st+subExtrapOpts.extrapFitWinPx;
    xFitX.Ind2 = length(posExtrapX.x)-subExtrapOpts.extrapPx1st...
        -subExtrapOpts.extrapFitWinPx+1:...
        length(posExtrapX.x)-subExtrapOpts.extrapPx1st;
    xFitX.xVals1 = posExtrapX.x(xFitX.Ind1);
    xFitX.xVals2 = posExtrapX.x(xFitX.Ind2);
    xFitX.extrap1 = posExtrapX.x(1:subExtrapOpts.extrapPx1st);
    xFitX.extrap2 = posExtrapX.x(end-subExtrapOpts.extrapPx1st+1:end);

    % Fit data structure for fitting X fields in the Y direction
    xFitY.Ind1 = 1+subExtrapOpts.extrapPx2nd:...
        subExtrapOpts.extrapPx2nd+subExtrapOpts.extrapFitWinPx;
    xFitY.Ind2 = length(posExtrapY.x)-subExtrapOpts.extrapPx2nd...
        -subExtrapOpts.extrapFitWinPx+1:...
        length(posExtrapY.x)-subExtrapOpts.extrapPx2nd;
    xFitY.xVals1 = posExtrapY.x(xFitY.Ind1);
    xFitY.xVals2 = posExtrapY.x(xFitY.Ind2);
    xFitY.extrap1 = posExtrapY.x(1:subExtrapOpts.extrapPx2nd);
    xFitY.extrap2 = posExtrapY.x(end-subExtrapOpts.extrapPx2nd+1:end);

    for tt = 1:st
        for pp = 1:ae_sy   
            % Create vectors of displacement values to fit, note: we have
            % already put the 'good' data into the extrap field
            xFitX.accelX1 = accelExtrap.x(pp,xFitX.Ind1,tt);
            xFitX.accelX2 = accelExtrap.x(pp,xFitX.Ind2,tt);
            if strcmp(subExtrapOpts.fieldOpt,'both')
                xFitY.dispY1 = accelExtrap.y(pp,xFitY.Ind1,tt);
                xFitY.dispY2 = accelExtrap.y(pp,xFitY.Ind2,tt);
            end

            if strcmp(subExtrapOpts.extrapMethod,'quadratic')
                % Extrapolate the 'X' edges (left and right hand sides of the
                % sample)
                accelExtrap.x(pp,1:subExtrapOpts.extrapPx1st,tt) = func_quick1DQuadraticExtrap(...
                    xFitX.xVals1,xFitX.accelX1,xFitX.extrap1);
                accelExtrap.x(pp,end-subExtrapOpts.extrapPx1st+1:end,tt) = func_quick1DQuadraticExtrap(...
                    xFitX.xVals2,xFitX.accelX2,xFitX.extrap2);

                % If we are cropping both edges then we need to extrapolate the
                % Y fields in the X direction as well
                if strcmp(subExtrapOpts.fieldOpt,'both')
                    accelExtrap.y(pp,1:subExtrapOpts.extrapPx2nd,tt) = func_quick1DQuadraticExtrap(...
                        xFitY.xVals1,xFitY.dispY1,xFitY.extrap1);
                    accelExtrap.y(pp,end-subExtrapOpts.extrapPx2nd+1:end,tt) = func_quick1DQuadraticExtrap(...
                        xFitY.xVals2,xFitY.dispY2,xFitY.extrap2);
                end
            elseif strcmp(subExtrapOpts.extrapMethod,'linear')
                % Extrapolate the 'X' edges (left and right hand sides of the
                % sample)
                accelExtrap.x(pp,1:subExtrapOpts.extrapPx1st,tt) = func_quick1DLinearExtrap(...
                    xFitX.xVals1,xFitX.accelX1,xFitX.xVals1(1),xFitX.accelX1(1),xFitX.extrap1);
                accelExtrap.x(pp,end-subExtrapOpts.extrapPx1st+1:end,tt) = func_quick1DLinearExtrap(...
                    xFitX.xVals2,xFitX.accelX2,xFitX.xVals2(end),xFitX.accelX2(end),xFitX.extrap2);
                
                % If we are cropping both edges then we need to extrapolate the
                % Y fields in the X direction as well
                if strcmp(subExtrapOpts.fieldOpt,'both')
                    accelExtrap.y(pp,1:subExtrapOpts.extrapPx2nd,tt) = func_quick1DLinearExtrap(...
                        xFitY.xVals1,xFitY.dispY1,xFitY.xVals1(1),xFitY.dispY1(1),xFitY.extrap1);
                    accelExtrap.y(pp,end-subExtrapOpts.extrapPx2nd+1:end,tt) = func_quick1DLinearExtrap(...
                        xFitY.xVals2,xFitY.dispY2,xFitY.xVals2(end),xFitY.dispY2(end),xFitY.extrap2);
                end
            end
        end
    end
    % Clear temp vars to make sure extrapolation in the Y direction doesn't use
    % these
    if debug
    figure; imagesc(flipud(accelExtrap.x(:,:,dbf))); title('accelExtrap.x');
    axis image; colorbar; colormap jet; set(gca,'YDir','normal');
    figure; imagesc(flipud(accelExtrap.y(:,:,dbf))); title('accelExtrap.y');
    axis image; colorbar; colormap jet; set(gca,'YDir','normal');

    end
    
    %--------------------------------------------------------------------------
    % Y DIRECTION
    yFitY.Ind1 = 1+subExtrapOpts.extrapPx1st:...
        subExtrapOpts.extrapPx1st+subExtrapOpts.extrapFitWinPx;
    yFitY.Ind2 = length(posExtrapY.y)-subExtrapOpts.extrapPx1st-subExtrapOpts.extrapFitWinPx+1:...
        length(posExtrapY.y)-subExtrapOpts.extrapPx1st;
    yFitY.yVals1 = posExtrapY.y(yFitY.Ind1);
    yFitY.yVals2 = posExtrapY.y(yFitY.Ind2);
    yFitY.extrap1 = posExtrapY.y(1:subExtrapOpts.extrapPx1st);
    yFitY.extrap2 = posExtrapY.y(end-subExtrapOpts.extrapPx1st+1:end);

    yFitX.Ind1 = 1+subExtrapOpts.extrapPx2nd:...
        subExtrapOpts.extrapPx2nd+subExtrapOpts.extrapFitWinPx;
    yFitX.Ind2 = length(posExtrapX.y)-subExtrapOpts.extrapPx2nd-subExtrapOpts.extrapFitWinPx+1:...
        length(posExtrapX.y)-subExtrapOpts.extrapPx2nd;
    yFitX.yVals1 = posExtrapX.y(yFitX.Ind1);
    yFitX.yVals2 = posExtrapX.y(yFitX.Ind2);
    yFitX.extrap1 = posExtrapX.y(1:subExtrapOpts.extrapPx2nd);
    yFitX.extrap2 = posExtrapX.y(end-subExtrapOpts.extrapPx2nd+1:end);

    for tt = 1:st
        for pp = 1:ae_sx
            % Create vectors of displacement values to fit, note: we have
            % already put the 'good' data into the extrap field
            yFitY.accelY1 = accelExtrap.y(yFitY.Ind1,pp,tt)';
            yFitY.accelY2 = accelExtrap.y(yFitY.Ind2,pp,tt)';
            if strcmp(subExtrapOpts.fieldOpt,'both')
                yFitX.accelX1 = accelExtrap.x(yFitX.Ind1,pp,tt)';
                yFitX.accelX2 = accelExtrap.x(yFitX.Ind2,pp,tt)';    
            end

            if strcmp(subExtrapOpts.extrapMethod,'quadratic')
                % Extrapolate the 'Y' edges (top and bottom sides of the
                % sample)
                accelExtrap.y(1:subExtrapOpts.extrapPx1st,pp,tt) = func_quick1DQuadraticExtrap(...
                    yFitY.yVals1,yFitY.accelY1,yFitY.extrap1);
                accelExtrap.y(end-subExtrapOpts.extrapPx1st+1:end,pp,tt) = func_quick1DQuadraticExtrap(...
                    yFitY.yVals2,yFitY.accelY2,yFitY.extrap2);

                % If we are cropping both edges then we need to extrapolate the
                % X fields in the Y direction as well
                if strcmp(subExtrapOpts.fieldOpt,'both')
                    accelExtrap.x(1:subExtrapOpts.extrapPx2nd,pp,tt) = func_quick1DQuadraticExtrap(...
                        yFitX.yVals1,yFitX.accelX1,yFitX.extrap1);
                    accelExtrap.x(end-subExtrapOpts.extrapPx2nd+1:end,pp,tt) = func_quick1DQuadraticExtrap(...
                        yFitX.yVals2,yFitX.accelX2,yFitX.extrap2);
                end
            elseif strcmp(subExtrapOpts.extrapMethod,'linear')
                % Extrapolate the 'Y' edges (top and bottom sides of the
                % sample)
                accelExtrap.y(1:subExtrapOpts.extrapPx1st,pp,tt) = func_quick1DLinearExtrap(...
                    yFitY.yVals1,yFitY.accelY1,yFitY.yVals1(1),yFitY.accelY1(1),yFitY.extrap1);
                accelExtrap.y(end-subExtrapOpts.extrapPx1st+1:end,pp,tt) = func_quick1DLinearExtrap(...
                    yFitY.yVals2,yFitY.accelY2,yFitY.yVals2(end),yFitY.accelY2(end),yFitY.extrap2);
                
                % If we are cropping both edges then we need to extrapolate the
                % X fields in the Y direction as well
                if strcmp(subExtrapOpts.fieldOpt,'both')
                    accelExtrap.x(1:subExtrapOpts.extrapPx2nd,pp,tt) = func_quick1DLinearExtrap(...
                        yFitX.yVals1,yFitX.accelX1,yFitX.yVals1(1),yFitX.accelX1(1),yFitX.extrap1);
                    accelExtrap.x(end-subExtrapOpts.extrapPx2nd+1:end,pp,tt) = func_quick1DLinearExtrap(...
                        yFitX.yVals2,yFitX.accelX2,yFitX.yVals2(end),yFitX.accelX2(end),yFitX.extrap2);
                end
            end
        end
    end

    if debug
    figure; imagesc(flipud(accelExtrap.x(:,:,dbf))); title('accelExtrap.x');
    axis image; colorbar; colormap jet; set(gca,'YDir','normal');
    figure; imagesc(flipud(accelExtrap.y(:,:,dbf))); title('accelExtrap.y');
    axis image; colorbar; colormap jet; set(gca,'YDir','normal');
    end

end

