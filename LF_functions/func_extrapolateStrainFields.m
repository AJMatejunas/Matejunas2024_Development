function strainExtrap = func_extrapolateStrainFields(...
    subExtrapOpts,pos,strainCrop,posCropX,posCropY,posCropS,debug)
% Extrapolates the missing edge data in the displacement fields for the
% grid method or DIC.
%
% Author: Lloyd Fletcher
% Image-Based Mechanics Group (IBMG), University of Southampton
% Date Created: 29/7/2020
% Date Edited: 31/7/2020
    
    % Turn off debug mode by default.
    if nargin < 6
        debug = false;
    end
    
    % Get the sizes of the fields for later calculations
    [ex_sy,ex_sx,st] = size(strainCrop.x);
    [ey_sy,ey_sx,~] = size(strainCrop.y);
    [es_sy,es_sx,~] = size(strainCrop.s);
    
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
        strainExtrap.x = nan(ex_sy,ex_sx+2*subExtrapOpts.extrapPx1st,st);
        strainExtrap.x(:,subExtrapOpts.extrapPx1st+1:end-subExtrapOpts.extrapPx1st,:) = strainCrop.x;
        strainExtrap.y = nan(ey_sy+2*subExtrapOpts.extrapPx1st,ey_sx,st);
        strainExtrap.y(subExtrapOpts.extrapPx1st+1:end-subExtrapOpts.extrapPx1st,:,:) = strainCrop.y;
        strainExtrap.s = nan(es_sy+2*subExtrapOpts.extrapPx1st,es_sx+2*subExtrapOpts.extrapPx1st,st);
        strainExtrap.s(subExtrapOpts.extrapPx1st+1:end-subExtrapOpts.extrapPx1st,...
            subExtrapOpts.extrapPx1st+1:end-subExtrapOpts.extrapPx1st,:) = strainCrop.s;
    else
        % Here we extrapolate in both directions, creates nans on all borders
        strainExtrap.x = nan(ex_sy+2*subExtrapOpts.extrapPx2nd,...
            ex_sx+2*subExtrapOpts.extrapPx1st,st);
        strainExtrap.x(1+subExtrapOpts.extrapPx2nd:end-subExtrapOpts.extrapPx2nd,...
            1+subExtrapOpts.extrapPx1st:end-subExtrapOpts.extrapPx1st,:) = strainCrop.x;

        strainExtrap.y = nan(ey_sy+2*subExtrapOpts.extrapPx1st,...
            ey_sx+2*subExtrapOpts.extrapPx2nd,st);
        strainExtrap.y(1+subExtrapOpts.extrapPx1st:end-subExtrapOpts.extrapPx1st,...
            1+subExtrapOpts.extrapPx2nd:end-subExtrapOpts.extrapPx2nd,:) = strainCrop.y;
        
        strainExtrap.s = nan(es_sy+2*subExtrapOpts.extrapPx1st,...
            es_sx+2*subExtrapOpts.extrapPx1st,st);
        strainExtrap.s(1+subExtrapOpts.extrapPx1st:end-subExtrapOpts.extrapPx1st,...
            1+subExtrapOpts.extrapPx1st:end-subExtrapOpts.extrapPx1st,:) = strainCrop.s;
    end
    [ee_sy,ee_sx,ee_st] = size(strainExtrap.x);
    
    if debug
    dbf = round(ee_st/2);    
    figure; imagesc(flipud(strainExtrap.x(:,:,dbf))); title('strainExtrap.x');
    axis image; colorbar; colormap jet; set(gca,'YDir','normal');
    figure; imagesc(flipud(strainExtrap.y(:,:,dbf))); title('strainExtrap.y');
    axis image; colorbar; colormap jet; set(gca,'YDir','normal');
    figure; imagesc(flipud(strainExtrap.s(:,:,dbf))); title('strainExtrap.s');
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
    
    % Extrapolate the same as in the primary direction for the shear.
    xFitS = xFitX;

    for tt = 1:st
        for pp = 1:ee_sy   
            % Create vectors of displacement values to fit, note: we have
            % already put the 'good' data into the extrap field
            xFitX.strainX1 = strainExtrap.x(pp,xFitX.Ind1,tt);
            xFitX.strainX2 = strainExtrap.x(pp,xFitX.Ind2,tt);
            xFitS.strainS1 = strainExtrap.s(pp,xFitS.Ind1,tt);
            xFitS.strainS2 = strainExtrap.s(pp,xFitS.Ind2,tt);
            if strcmp(subExtrapOpts.fieldOpt,'both')
                xFitY.dispY1 = strainExtrap.y(pp,xFitY.Ind1,tt);
                xFitY.dispY2 = strainExtrap.y(pp,xFitY.Ind2,tt);
            end

            if strcmp(subExtrapOpts.extrapMethod,'quadratic')
                % Extrapolate the 'X' edges (left and right hand sides of the
                % sample)
                strainExtrap.x(pp,1:subExtrapOpts.extrapPx1st,tt) = func_quick1DQuadraticExtrap(...
                    xFitX.xVals1,xFitX.strainX1,xFitX.extrap1);
                strainExtrap.x(pp,end-subExtrapOpts.extrapPx1st+1:end,tt) = func_quick1DQuadraticExtrap(...
                    xFitX.xVals2,xFitX.strainX2,xFitX.extrap2);
                
                % Shear strain treated the same as the primary direction
                strainExtrap.s(pp,1:subExtrapOpts.extrapPx1st,tt) = func_quick1DQuadraticExtrap(...
                    xFitS.xVals1,xFitS.strainS1,xFitS.extrap1);
                strainExtrap.s(pp,end-subExtrapOpts.extrapPx1st+1:end,tt) = func_quick1DQuadraticExtrap(...
                    xFitS.xVals2,xFitS.strainS2,xFitS.extrap2);

                % If we are cropping both edges then we need to extrapolate the
                % Y fields in the X direction as well
                if strcmp(subExtrapOpts.fieldOpt,'both')
                    strainExtrap.y(pp,1:subExtrapOpts.extrapPx2nd,tt) = func_quick1DQuadraticExtrap(...
                        xFitY.xVals1,xFitY.dispY1,xFitY.extrap1);
                    strainExtrap.y(pp,end-subExtrapOpts.extrapPx2nd+1:end,tt) = func_quick1DQuadraticExtrap(...
                        xFitY.xVals2,xFitY.dispY2,xFitY.extrap2);
                end
            elseif strcmp(subExtrapOpts.extrapMethod,'linear')
                % Extrapolate the 'X' edges (left and right hand sides of the
                % sample)
                strainExtrap.x(pp,1:subExtrapOpts.extrapPx1st,tt) = func_quick1DLinearExtrap(...
                    xFitX.xVals1,xFitX.strainX1,xFitX.xVals1(1),xFitX.strainX1(1),xFitX.extrap1);
                strainExtrap.x(pp,end-subExtrapOpts.extrapPx1st+1:end,tt) = func_quick1DLinearExtrap(...
                    xFitX.xVals2,xFitX.strainX2,xFitX.xVals2(end),xFitX.strainX2(end),xFitX.extrap2);
                
                % Shear strain treated the same as the primary direction
                strainExtrap.s(pp,1:subExtrapOpts.extrapPx1st,tt) = func_quick1DLinearExtrap(...
                    xFitS.xVals1,xFitS.strainS1,xFitS.xVals1(1),xFitS.strainS1(1),xFitS.extrap1);
                strainExtrap.s(pp,end-subExtrapOpts.extrapPx1st+1:end,tt) = func_quick1DLinearExtrap(...
                    xFitS.xVals2,xFitS.strainS2,xFitS.xVals2(end),xFitS.strainS2(end),xFitS.extrap2);

                % If we are cropping both edges then we need to extrapolate the
                % Y fields in the X direction as well
                if strcmp(subExtrapOpts.fieldOpt,'both')
                    strainExtrap.y(pp,1:subExtrapOpts.extrapPx2nd,tt) = func_quick1DLinearExtrap(...
                        xFitY.xVals1,xFitY.dispY1,xFitY.xVals1(1),xFitY.dispY1(1),xFitY.extrap1);
                    strainExtrap.y(pp,end-subExtrapOpts.extrapPx2nd+1:end,tt) = func_quick1DLinearExtrap(...
                        xFitY.xVals2,xFitY.dispY2,xFitY.xVals2(end),xFitY.dispY2(end),xFitY.extrap2);
                end
            end
        end
    end
    % Clear temp vars to make sure extrapolation in the Y direction doesn't use
    % these
    if debug
    figure; imagesc(flipud(strainExtrap.x(:,:,dbf))); title('strainExtrap.x');
    axis image; colorbar; colormap jet; set(gca,'YDir','normal');
    figure; imagesc(flipud(strainExtrap.y(:,:,dbf))); title('strainExtrap.y');
    axis image; colorbar; colormap jet; set(gca,'YDir','normal');
    figure; imagesc(flipud(strainExtrap.s(:,:,dbf))); title('strainExtrap.s');
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
    
    % Extrapolate the same as in the primary direction for the shear.
    yFitS = yFitY;

    for tt = 1:st
        for pp = 1:ee_sx
            % Create vectors of displacement values to fit, note: we have
            % already put the 'good' data into the extrap field
            yFitY.strainY1 = strainExtrap.y(yFitY.Ind1,pp,tt)';
            yFitY.strainY2 = strainExtrap.y(yFitY.Ind2,pp,tt)';
            % Treat the shear the same as the primary direction
            yFitS.strainS1 = strainExtrap.s(yFitS.Ind1,pp,tt)';
            yFitS.strainS2 = strainExtrap.s(yFitS.Ind2,pp,tt)';
            if strcmp(subExtrapOpts.fieldOpt,'both')
                yFitX.strainX1 = strainExtrap.x(yFitX.Ind1,pp,tt)';
                yFitX.strainX2 = strainExtrap.x(yFitX.Ind2,pp,tt)';    
            end

            if strcmp(subExtrapOpts.extrapMethod,'quadratic')
                % Extrapolate the 'Y' edges (top and bottom sides of the
                % sample)
                strainExtrap.y(1:subExtrapOpts.extrapPx1st,pp,tt) = func_quick1DQuadraticExtrap(...
                    yFitY.yVals1,yFitY.strainY1,yFitY.extrap1);
                strainExtrap.y(end-subExtrapOpts.extrapPx1st+1:end,pp,tt) = func_quick1DQuadraticExtrap(...
                    yFitY.yVals2,yFitY.strainY2,yFitY.extrap2);
                
                % Treat the shear the same as the primary direction
                strainExtrap.s(1:subExtrapOpts.extrapPx1st,pp,tt) = func_quick1DQuadraticExtrap(...
                    yFitS.yVals1,yFitS.strainS1,yFitS.extrap1);
                strainExtrap.s(end-subExtrapOpts.extrapPx1st+1:end,pp,tt) = func_quick1DQuadraticExtrap(...
                    yFitS.yVals2,yFitS.strainS2,yFitS.extrap2);

                % If we are cropping both edges then we need to extrapolate the
                % X fields in the Y direction as well
                if strcmp(subExtrapOpts.fieldOpt,'both')
                    strainExtrap.x(1:subExtrapOpts.extrapPx2nd,pp,tt) = func_quick1DQuadraticExtrap(...
                        yFitX.yVals1,yFitX.strainX1,yFitX.extrap1);
                    strainExtrap.x(end-subExtrapOpts.extrapPx2nd+1:end,pp,tt) = func_quick1DQuadraticExtrap(...
                        yFitX.yVals2,yFitX.strainX2,yFitX.extrap2);
                end
            elseif strcmp(subExtrapOpts.extrapMethod,'linear')
                % Extrapolate the 'Y' edges (top and bottom sides of the
                % sample)
                strainExtrap.y(1:subExtrapOpts.extrapPx1st,pp,tt) = func_quick1DLinearExtrap(...
                    yFitY.yVals1,yFitY.strainY1,yFitY.yVals1(1),yFitY.strainY1(1),yFitY.extrap1);
                strainExtrap.y(end-subExtrapOpts.extrapPx1st+1:end,pp,tt) = func_quick1DLinearExtrap(...
                    yFitY.yVals2,yFitY.strainY2,yFitY.yVals2(end),yFitY.strainY2(end),yFitY.extrap2);
                
                % Treat the shear the same as the primary direction
                strainExtrap.s(1:subExtrapOpts.extrapPx1st,pp,tt) = func_quick1DLinearExtrap(...
                    yFitS.yVals1,yFitS.strainS1,yFitS.yVals1(1),yFitS.strainS1(1),yFitS.extrap1);
                strainExtrap.s(end-subExtrapOpts.extrapPx1st+1:end,pp,tt) = func_quick1DLinearExtrap(...
                    yFitS.yVals2,yFitS.strainS2,yFitS.yVals2(end),yFitS.strainS2(end),yFitS.extrap2);
                
                % If we are cropping both edges then we need to extrapolate the
                % X fields in the Y direction as well
                if strcmp(subExtrapOpts.fieldOpt,'both')
                    strainExtrap.x(1:subExtrapOpts.extrapPx2nd,pp,tt) = func_quick1DLinearExtrap(...
                        yFitX.yVals1,yFitX.strainX1,yFitX.yVals1(1),yFitX.strainX1(1),yFitX.extrap1);
                    strainExtrap.x(end-subExtrapOpts.extrapPx2nd+1:end,pp,tt) = func_quick1DLinearExtrap(...
                        yFitX.yVals2,yFitX.strainX2,yFitX.yVals2(end),yFitX.strainX2(end),yFitX.extrap2);
                end
            end
        end
    end

    if debug
    figure; imagesc(flipud(strainExtrap.x(:,:,dbf))); title('strainExtrap.x');
    axis image; colorbar; colormap jet; set(gca,'YDir','normal');
    figure; imagesc(flipud(strainExtrap.y(:,:,dbf))); title('strainExtrap.y');
    axis image; colorbar; colormap jet; set(gca,'YDir','normal');
    figure; imagesc(flipud(strainExtrap.s(:,:,dbf))); title('strainExtrap.s');
    axis image; colorbar; colormap jet; set(gca,'YDir','normal');
    end

end

