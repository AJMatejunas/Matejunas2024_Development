function [extrapVar] = func_spatialExtrapFullFieldVar(pos,inputVar,extrapDim,...
    extrapPixels,extrapMethod)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 27/3/2017
% Extrapolate full field variable along the specified dimension
% Defaults to extrapolating along the first dimension

    % Delete the edge pixels and assign the position variable based on the
    % dimension we are extrapolating
    if extrapDim == 2
        stepDim = 1;
        crop_pos = pos.x(extrapPixels+1:end-extrapPixels);
        crop_inputVar = inputVar(:,extrapPixels+1:end-extrapPixels,:);
    else
        stepDim = 2;
        crop_pos = pos.y(extrapPixels+1:end-extrapPixels);
        crop_inputVar = inputVar(extrapPixels+1:end-extrapPixels,:,:);     
    end
    
    % Fill in everything in the variable to be extrapolated, fill in the
    % missing pitch later
    extrapVar = inputVar;
    temp_extrapVar = zeros(size(inputVar));
    
    % If the option is nearest then we just need to pad the edges of the
    % data
    if strcmp(extrapMethod,'nearest')
        if extrapDim == 2
            for t = 1:size(inputVar,3)
                for p = 1:size(inputVar,stepDim)
                    % Add the extrapolated data back onto the edges
                    extrapVar(p,1:extrapPixels,t) = crop_inputVar(p,1,t);
                    extrapVar(p,end-extrapPixels+1:end,t) = crop_inputVar(p,end,t);
                end
            end
        else
            for t = 1:size(inputVar,3)
                for p = 1:size(inputVar,stepDim)
                    % Pad the data with the edge values
                    extrapVar(1:extrapPixels,p,t) = crop_inputVar(1,p,t);
                    extrapVar(end-extrapPixels+1:end,p,t) = crop_inputVar(end,p,t);
                end
            end     
        end
    else
        if extrapDim == 2
            for t = 1:size(inputVar,3)
                for p = 1:size(inputVar,stepDim)
                    % Extrapolate the variable over the given number of pixels
                    temp_extrapVar(p,:,t) = interp1(crop_pos,crop_inputVar(p,:,t),pos.x,extrapMethod,'extrap');
                    % Add the extrapolated data back onto the edges
                    extrapVar(p,1:extrapPixels,t) = temp_extrapVar(p,1:extrapPixels,t);
                    extrapVar(p,end-extrapPixels+1:end,t) = temp_extrapVar(p,end-extrapPixels+1:end,t);
                end
            end
        else
            for t = 1:size(inputVar,3)
                for p = 1:size(inputVar,stepDim)
                    % Extrapolate the variable over the given number of pixels
                    temp_extrapVar(:,p,t) = interp1(crop_pos,crop_inputVar(:,p,t),pos.y,extrapMethod,'extrap');
                    % Add the extrapolated data back onto the edges
                    extrapVar(1:extrapPixels,p,t) = temp_extrapVar(1:extrapPixels,p,t);
                    extrapVar(end-extrapPixels+1:end,p,t) = temp_extrapVar(end-extrapPixels+1:end,p,t);
                end
            end     
        end
    end
end

