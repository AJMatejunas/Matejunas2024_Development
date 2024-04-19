function [disp,vel,accel,strain] = func_spatialExtrapAvgDataOpts(extrapOpts,time,pos,disp,vel,accel,strain)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 9/2/2017
%
% Spatial extrapolation of average data

    % Delete the edge pixels
    temp_posX = pos.x(extrapOpts.dispPx+1:end-extrapOpts.dispPx);
    temp_dispXAvg = disp.xAvg(extrapOpts.dispPx+1:end-extrapOpts.dispPx,:);
    temp_velXAvg = vel.xAvg(extrapOpts.dispPx+1:end-extrapOpts.dispPx,:);
    temp_accelXAvg = accel.xAvg(extrapOpts.dispPx+1:end-extrapOpts.dispPx,:); 
    
    temp_strainPosX = pos.x(extrapOpts.strainPx+1:end-extrapOpts.strainPx);
    temp_strainXAvg = strain.xAvg(extrapOpts.strainPx+1:end-extrapOpts.strainPx,:);
    temp_strainYAvg = strain.yAvg(extrapOpts.strainPx+1:end-extrapOpts.strainPx,:);
     
    for i = 1:time.numFrames
        % Extrapolate displacement, velocity and acceleration
        temp_dispXAvg_ext(:,i) = interp1(temp_posX,temp_dispXAvg(:,i),pos.x,extrapOpts.dispMethod,'extrap');
        disp.xAvg(1:extrapOpts.dispPx,i) = temp_dispXAvg_ext(1:extrapOpts.dispPx,i);
        disp.xAvg(end-extrapOpts.dispPx+1:end,i) = temp_dispXAvg_ext(end-extrapOpts.dispPx+1:end,i);
        
        temp_velXAvg_ext(:,i) = interp1(temp_posX,temp_velXAvg(:,i),pos.x,extrapOpts.dispMethod,'extrap');
        vel.xAvg(1:extrapOpts.dispPx,i) = temp_velXAvg_ext(1:extrapOpts.dispPx,i);
        vel.xAvg(end-extrapOpts.dispPx+1:end,i) = temp_velXAvg_ext(end-extrapOpts.dispPx+1:end,i);
        
        temp_accelXAvg_ext(:,i) = interp1(temp_posX,temp_accelXAvg(:,i),pos.x,extrapOpts.dispMethod,'extrap');
        accel.xAvg(1:extrapOpts.dispPx,i) = temp_accelXAvg_ext(1:extrapOpts.dispPx,i);
        accel.xAvg(end-extrapOpts.dispPx+1:end,i) = temp_accelXAvg_ext(end-extrapOpts.dispPx+1:end,i);
        
        % Extrapolate the strains
        temp_strainXAvg_ext(:,i) = interp1(temp_strainPosX,temp_strainXAvg(:,i),pos.x,extrapOpts.strainMethod,'extrap');
        strain.xAvg(1:extrapOpts.strainPx,i) = temp_strainXAvg_ext(1:extrapOpts.strainPx,i);
        strain.xAvg(end-extrapOpts.strainPx+1:end,i) = temp_strainXAvg_ext(end-extrapOpts.strainPx+1:end,i);
        
        temp_strainYAvg_ext(:,i) = interp1(temp_strainPosX,temp_strainYAvg(:,i),pos.x,extrapOpts.strainMethod,'extrap');
        strain.yAvg(1:extrapOpts.strainPx,i) = temp_strainYAvg_ext(1:extrapOpts.strainPx,i);
        strain.yAvg(end-extrapOpts.strainPx+1:end,i) = temp_strainYAvg_ext(end-extrapOpts.strainPx+1:end,i);
    end

end

