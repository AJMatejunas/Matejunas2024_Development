function [disp,vel,accel,strain] = func_spatialExtrapAvgData(time,pos,disp,vel,accel,strain,method,dispExtPx,strainExtPx)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 9/2/2017
%
% Spatial extrapolation of average data

    % Delete the edge pixels
    temp_posX = pos.x(dispExtPx+1:end-dispExtPx);
    temp_dispXAvg = disp.xAvg(dispExtPx+1:end-dispExtPx,:);
    temp_velXAvg = vel.xAvg(dispExtPx+1:end-dispExtPx,:);
    temp_accelXAvg = accel.xAvg(dispExtPx+1:end-dispExtPx,:); 
    
    temp_strainPosX = pos.x(strainExtPx+1:end-strainExtPx);
    temp_strainXAvg = strain.xAvg(strainExtPx+1:end-strainExtPx,:);
    temp_strainYAvg = strain.yAvg(strainExtPx+1:end-strainExtPx,:);
     
    for i = 1:time.numFrames
        % Extrapolate displacement, velocity and acceleration
        temp_dispXAvg_ext(:,i) = interp1(temp_posX,temp_dispXAvg(:,i),pos.x,method,'extrap');
        disp.xAvg(1:dispExtPx,i) = temp_dispXAvg_ext(1:dispExtPx,i);
        disp.xAvg(end-dispExtPx+1:end,i) = temp_dispXAvg_ext(end-dispExtPx+1:end,i);
        
        temp_velXAvg_ext(:,i) = interp1(temp_posX,temp_velXAvg(:,i),pos.x,method,'extrap');
        vel.xAvg(1:dispExtPx,i) = temp_velXAvg_ext(1:dispExtPx,i);
        vel.xAvg(end-dispExtPx+1:end,i) = temp_velXAvg_ext(end-dispExtPx+1:end,i);
        
        temp_accelXAvg_ext(:,i) = interp1(temp_posX,temp_accelXAvg(:,i),pos.x,method,'extrap');
        accel.xAvg(1:dispExtPx,i) = temp_accelXAvg_ext(1:dispExtPx,i);
        accel.xAvg(end-dispExtPx+1:end,i) = temp_accelXAvg_ext(end-dispExtPx+1:end,i);
        
        % Extrapolate the strains
        temp_strainXAvg_ext(:,i) = interp1(temp_strainPosX,temp_strainXAvg(:,i),pos.x,method,'extrap');
        strain.xAvg(1:strainExtPx,i) = temp_strainXAvg_ext(1:strainExtPx,i);
        strain.xAvg(end-strainExtPx+1:end,i) = temp_strainXAvg_ext(end-strainExtPx+1:end,i);
        
        temp_strainYAvg_ext(:,i) = interp1(temp_strainPosX,temp_strainYAvg(:,i),pos.x,method,'extrap');
        strain.yAvg(1:strainExtPx,i) = temp_strainYAvg_ext(1:strainExtPx,i);
        strain.yAvg(end-strainExtPx+1:end,i) = temp_strainYAvg_ext(end-strainExtPx+1:end,i);
    end

end

