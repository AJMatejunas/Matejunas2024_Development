function disp = func_padWidthPixels(pos,disp,padWidthPixels,padMethod)
% Pad over the width to remove nan values from method 2 
    % Delete the edge pixels
    temp_posY = pos.y(padWidthPixels+1:end-padWidthPixels);
    temp_dispX = disp.x(padWidthPixels+1:end-padWidthPixels,:,:);   
    
    % Use a linear extrapolation or padded extrapolation
    for f = 1:size(disp.x,3)
        for p = 1:size(disp.x,2)
            % Extrapolate displacement, velocity and acceleration
            temp_dispX_ext(:,p,f) = interp1(temp_posY,temp_dispX(:,p,f),pos.y,padMethod,'extrap');
            disp.x(1:padWidthPixels,p,f) = temp_dispX_ext(1:padWidthPixels,p,f);
            disp.x(end-padWidthPixels+1:end,p,f) = temp_dispX_ext(end-padWidthPixels+1:end,p,f);
        end
    end

end

