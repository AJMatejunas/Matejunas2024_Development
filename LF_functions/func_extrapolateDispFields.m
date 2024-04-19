function dispExtrap = func_extrapolateDispFields(...
    extrapOptsD,pos,dispCrop,posCropX,posCropY,debug)
% Extrapolates the missing edge data in the displacement fields for the
% grid method or DIC.
%
% Author: Lloyd Fletcher
% Image-Based Mechanics Group (IBMG), University of Southampton
% Date Created: 28/7/2020

    % Turn off debug mode by default.
    if nargin < 6
        debug = false;
    end

    % Get the sizes of the fields for later calculations
    [ux_sy,ux_sx,st] = size(dispCrop.x);
    [uy_sy,uy_sx,~] = size(dispCrop.y);

    % Define extrapolation co-ordinates - X displacement field
    posExtrapX.x = min(posCropX.x)-extrapOptsD.extrapPx1st*pos.xStep:pos.xStep:...
        max(posCropX.x)+extrapOptsD.extrapPx1st*pos.xStep;
    posExtrapX.y = min(posCropX.y)-extrapOptsD.extrapPx2nd*pos.yStep:pos.yStep:...
        max(posCropX.y)+extrapOptsD.extrapPx2nd*pos.yStep;
    % Define extrapolation co-ordinates - Y displacement field
    posExtrapY.x = min(posCropY.x)-extrapOptsD.extrapPx2nd*pos.xStep:pos.xStep:...
        max(posCropY.x)+extrapOptsD.extrapPx2nd*pos.xStep;
    posExtrapY.y = min(posCropY.y)-extrapOptsD.extrapPx1st*pos.yStep:pos.yStep:...
        max(posCropY.y)+extrapOptsD.extrapPx1st*pos.yStep;

    % Initiliase 'dispExtrap' vars for speed
    if strcmp(extrapOptsD.fieldOpt,'XY')
        % Here we only extrapolate X in the X direction and same for Y
        dispExtrap.x = nan(ux_sy,ux_sx+2*extrapOptsD.extrapPx1st,st);
        dispExtrap.x(:,extrapOptsD.extrapPx1st+1:end-extrapOptsD.extrapPx1st,:) = dispCrop.x;
        dispExtrap.y = nan(uy_sy+2*extrapOptsD.extrapPx1st,uy_sx,st);
        dispExtrap.y(extrapOptsD.extrapPx1st+1:end-extrapOptsD.extrapPx1st,:,:) = dispCrop.y;
    else
        % Here we extrapolate in both directions, creates zeros on all borders
        dispExtrap.x = nan(ux_sy+2*extrapOptsD.extrapPx2nd,...
            ux_sx+2*extrapOptsD.extrapPx1st,st);
        dispExtrap.x(1+extrapOptsD.extrapPx2nd:end-extrapOptsD.extrapPx2nd,...
            1+extrapOptsD.extrapPx1st:end-extrapOptsD.extrapPx1st,:) = dispCrop.x;

        dispExtrap.y = nan(uy_sy+2*extrapOptsD.extrapPx1st,...
            uy_sx+2*extrapOptsD.extrapPx2nd,st);
        dispExtrap.y(1+extrapOptsD.extrapPx1st:end-extrapOptsD.extrapPx1st,...
            1+extrapOptsD.extrapPx2nd:end-extrapOptsD.extrapPx2nd,:) = dispCrop.y;
    end
    [ue_sy,ue_sx,~] = size(dispExtrap.x);

    if debug
    figure; imagesc(flipud(dispExtrap.x(:,:,end))); title('dispExtrap.x');
    axis image; colorbar; colormap jet; set(gca,'YDir','normal');
    figure; imagesc(flipud(dispExtrap.y(:,:,end))); title('dispExtrap.y');
    axis image; colorbar; colormap jet; set(gca,'YDir','normal');
    end

    %--------------------------------------------------------------------------
    % X DIRECTION
    % Fit data structure for fitting X fields in the X direction
    xFitX.Ind1 = 1+extrapOptsD.extrapPx1st:...
        extrapOptsD.extrapPx1st+extrapOptsD.extrapFitWinPx;
    xFitX.Ind2 = length(posExtrapX.x)-extrapOptsD.extrapPx1st...
        -extrapOptsD.extrapFitWinPx+1:...
        length(posExtrapX.x)-extrapOptsD.extrapPx1st;
    xFitX.xVals1 = posExtrapX.x(xFitX.Ind1);
    xFitX.xVals2 = posExtrapX.x(xFitX.Ind2);
    xFitX.extrap1 = posExtrapX.x(1:extrapOptsD.extrapPx1st);
    xFitX.extrap2 = posExtrapX.x(end-extrapOptsD.extrapPx1st+1:end);

    % Fit data structure for fitting X fields in the Y direction
    xFitY.Ind1 = 1+extrapOptsD.extrapPx2nd:...
        extrapOptsD.extrapPx2nd+extrapOptsD.extrapFitWinPx;
    xFitY.Ind2 = length(posExtrapY.x)-extrapOptsD.extrapPx2nd...
        -extrapOptsD.extrapFitWinPx+1:...
        length(posExtrapY.x)-extrapOptsD.extrapPx2nd;
    xFitY.xVals1 = posExtrapY.x(xFitY.Ind1);
    xFitY.xVals2 = posExtrapY.x(xFitY.Ind2);
    xFitY.extrap1 = posExtrapY.x(1:extrapOptsD.extrapPx2nd);
    xFitY.extrap2 = posExtrapY.x(end-extrapOptsD.extrapPx2nd+1:end);

    for tt = 1:st
        for pp = 1:ue_sy   
            % Create vectors of displacement values to fit, note: we have
            % already put the 'good' data into the extrap field
            xFitX.dispX1 = dispExtrap.x(pp,xFitX.Ind1,tt);
            xFitX.dispX2 = dispExtrap.x(pp,xFitX.Ind2,tt);
            if strcmp(extrapOptsD.fieldOpt,'both')
                xFitY.dispY1 = dispExtrap.y(pp,xFitY.Ind1,tt);
                xFitY.dispY2 = dispExtrap.y(pp,xFitY.Ind2,tt);
            end

            if strcmp(extrapOptsD.extrapMethod,'quadratic')
                % Extrapolate the 'X' edges (left and right hand sides of the
                % sample)
                dispExtrap.x(pp,1:extrapOptsD.extrapPx1st,tt) = func_quick1DQuadraticExtrap(...
                    xFitX.xVals1,xFitX.dispX1,xFitX.extrap1);
                dispExtrap.x(pp,end-extrapOptsD.extrapPx1st+1:end,tt) = func_quick1DQuadraticExtrap(...
                    xFitX.xVals2,xFitX.dispX2,xFitX.extrap2);

                % If we are cropping both edges then we need to extrapolate the
                % Y fields in the X direction as well
                if strcmp(extrapOptsD.fieldOpt,'both')
                    dispExtrap.y(pp,1:extrapOptsD.extrapPx2nd,tt) = func_quick1DQuadraticExtrap(...
                        xFitY.xVals1,xFitY.dispY1,xFitY.extrap1);
                    dispExtrap.y(pp,end-extrapOptsD.extrapPx2nd+1:end,tt) = func_quick1DQuadraticExtrap(...
                        xFitY.xVals2,xFitY.dispY2,xFitY.extrap2);
                end
            elseif strcmp(extrapOptsD.extrapMethod,'linear')
                % Extrapolate the 'X' edges (left and right hand sides of the
                % sample)
                dispExtrap.x(pp,1:extrapOptsD.extrapPx1st,tt) = func_quick1DLinearExtrap(...
                    xFitX.xVals1,xFitX.dispX1,xFitX.xVals1(1),xFitX.dispX1(1),xFitX.extrap1);
                dispExtrap.x(pp,end-extrapOptsD.extrapPx1st+1:end,tt) = func_quick1DLinearExtrap(...
                    xFitX.xVals2,xFitX.dispX2,xFitX.xVals2(end),xFitX.dispX2(end),xFitX.extrap2);

                % If we are cropping both edges then we need to extrapolate the
                % Y fields in the X direction as well
                if strcmp(extrapOptsD.fieldOpt,'both')
                    dispExtrap.y(pp,1:extrapOptsD.extrapPx2nd,tt) = func_quick1DLinearExtrap(...
                        xFitY.xVals1,xFitY.dispY1,xFitY.xVals1(1),xFitY.dispY1(1),xFitY.extrap1);
                    dispExtrap.y(pp,end-extrapOptsD.extrapPx2nd+1:end,tt) = func_quick1DLinearExtrap(...
                        xFitY.xVals2,xFitY.dispY2,xFitY.xVals2(end),xFitY.dispY2(end),xFitY.extrap2);
                end
            end
        end
    end
    % Clear temp vars to make sure extrapolation in the Y direction doesn't use
    % these
    if debug
    figure; imagesc(flipud(dispExtrap.x(:,:,end))); title('dispExtrap.x');
    axis image; colorbar; colormap jet; set(gca,'YDir','normal');
    figure; imagesc(flipud(dispExtrap.y(:,:,end))); title('dispExtrap.y');
    axis image; colorbar; colormap jet; set(gca,'YDir','normal');
    end

    %--------------------------------------------------------------------------
    % Y DIRECTION
    yFitY.Ind1 = 1+extrapOptsD.extrapPx1st:...
        extrapOptsD.extrapPx1st+extrapOptsD.extrapFitWinPx;
    yFitY.Ind2 = length(posExtrapY.y)-extrapOptsD.extrapPx1st-extrapOptsD.extrapFitWinPx+1:...
        length(posExtrapY.y)-extrapOptsD.extrapPx1st;
    yFitY.yVals1 = posExtrapY.y(yFitY.Ind1);
    yFitY.yVals2 = posExtrapY.y(yFitY.Ind2);
    yFitY.extrap1 = posExtrapY.y(1:extrapOptsD.extrapPx1st);
    yFitY.extrap2 = posExtrapY.y(end-extrapOptsD.extrapPx1st+1:end);

    yFitX.Ind1 = 1+extrapOptsD.extrapPx2nd:...
        extrapOptsD.extrapPx2nd+extrapOptsD.extrapFitWinPx;
    yFitX.Ind2 = length(posExtrapX.y)-extrapOptsD.extrapPx2nd-extrapOptsD.extrapFitWinPx+1:...
        length(posExtrapX.y)-extrapOptsD.extrapPx2nd;
    yFitX.yVals1 = posExtrapX.y(yFitX.Ind1);
    yFitX.yVals2 = posExtrapX.y(yFitX.Ind2);
    yFitX.extrap1 = posExtrapX.y(1:extrapOptsD.extrapPx2nd);
    yFitX.extrap2 = posExtrapX.y(end-extrapOptsD.extrapPx2nd+1:end);

    for tt = 1:st
        for pp = 1:ue_sx
            % Create vectors of displacement values to fit, note: we have
            % already put the 'good' data into the extrap field
            yFitY.dispY1 = dispExtrap.y(yFitY.Ind1,pp,tt)';
            yFitY.dispY2 = dispExtrap.y(yFitY.Ind2,pp,tt)';
            if strcmp(extrapOptsD.fieldOpt,'both')
                yFitX.dispX1 = dispExtrap.x(yFitX.Ind1,pp,tt)';
                yFitX.dispX2 = dispExtrap.x(yFitX.Ind2,pp,tt)';    
            end

            if strcmp(extrapOptsD.extrapMethod,'quadratic')
                % Extrapolate the 'Y' edges (top and bottom sides of the
                % sample)
                dispExtrap.y(1:extrapOptsD.extrapPx1st,pp,tt) = func_quick1DQuadraticExtrap(...
                    yFitY.yVals1,yFitY.dispY1,yFitY.extrap1);
                dispExtrap.y(end-extrapOptsD.extrapPx1st+1:end,pp,tt) = func_quick1DQuadraticExtrap(...
                    yFitY.yVals2,yFitY.dispY2,yFitY.extrap2);

                % If we are cropping both edges then we need to extrapolate the
                % X fields in the Y direction as well
                if strcmp(extrapOptsD.fieldOpt,'both')
                    dispExtrap.x(1:extrapOptsD.extrapPx2nd,pp,tt) = func_quick1DQuadraticExtrap(...
                        yFitX.yVals1,yFitX.dispX1,yFitX.extrap1);
                    dispExtrap.x(end-extrapOptsD.extrapPx2nd+1:end,pp,tt) = func_quick1DQuadraticExtrap(...
                        yFitX.yVals2,yFitX.dispX2,yFitX.extrap2);
                end
            elseif strcmp(extrapOptsD.extrapMethod,'linear')
                % Extrapolate the 'Y' edges (top and bottom sides of the
                % sample)
                dispExtrap.y(1:extrapOptsD.extrapPx1st,pp,tt) = func_quick1DLinearExtrap(...
                    yFitY.yVals1,yFitY.dispY1,yFitY.yVals1(1),yFitY.dispY1(1),yFitY.extrap1);
                dispExtrap.y(end-extrapOptsD.extrapPx1st+1:end,pp,tt) = func_quick1DLinearExtrap(...
                    yFitY.yVals2,yFitY.dispY2,yFitY.yVals2(end),yFitY.dispY2(end),yFitY.extrap2);

                % If we are cropping both edges then we need to extrapolate the
                % X fields in the Y direction as well
                if strcmp(extrapOptsD.fieldOpt,'both')
                    dispExtrap.x(1:extrapOptsD.extrapPx2nd,pp,tt) = func_quick1DLinearExtrap(...
                        yFitX.yVals1,yFitX.dispX1,yFitX.yVals1(1),yFitX.dispX1(1),yFitX.extrap1);
                    dispExtrap.x(end-extrapOptsD.extrapPx2nd+1:end,pp,tt) = func_quick1DLinearExtrap(...
                        yFitX.yVals2,yFitX.dispX2,yFitX.yVals2(end),yFitX.dispX2(end),yFitX.extrap2);
                end
            end
        end
    end

    if debug
    figure; imagesc(flipud(dispExtrap.x(:,:,end))); title('dispExtrap.x');
    axis image; colorbar; colormap jet; set(gca,'YDir','normal');
    figure; imagesc(flipud(dispExtrap.y(:,:,end))); title('dispExtrap.y');
    axis image; colorbar; colormap jet; set(gca,'YDir','normal');
    end
end

