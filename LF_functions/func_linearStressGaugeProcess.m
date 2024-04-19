function stress = func_linearStressGaugeProcess(specimen,material,pos,accel,stress)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 31/3/2017
% Calculates the first moment of the stress gauge and caculates the axial
% stress distrbution using a linear approximation

% Get the iteration limits
[sizeY,sizeX,sizeT] = size(accel.x);

% Set the co-ord system central to the specimen
posVecX = pos.x;
posVecY = pos.y - specimen.height/2;
[posXGrid,posYGrid] = meshgrid(posVecX,posVecY);
posXGridF = padarray(posXGrid,[0,0,sizeT-1],'replicate','post'); 
posYGridF = padarray(posYGrid,[0,0,sizeT-1],'replicate','post'); 

% Calculate the weighted averages
ax_y = accel.x.*posYGridF;
ay_x = accel.y.*posXGridF;

% Pre-alloc for speed
avg_ax_y = zeros(sizeX,sizeT);
avg_ay_x = zeros(sizeX,sizeT);
avg_ay = zeros(sizeX,sizeT);

    
for x = 1:sizeX
    % Calculate the surface average acceleration up to 'x'
    avg_ax_y(x,:) = mean(mean(ax_y(:,1:x,:)));
    avg_ay_x(x,:) = mean(mean(ay_x(:,1:x,:)));
    avg_ay(x,:) = mean(mean(accel.y(:,1:x,:)));

    % Calculate the first moment and the linear distribution of stress
    stress.xFirstMoment(x,:) = (12*material.rho*posVecX(x)/specimen.height)*...
        (avg_ax_y(x,:)-avg_ay_x(x,:)+posVecX(x)*avg_ay(x,:));

    for y = 1:sizeY
        stress.xLinearGauge(y,x,:) = stress.xAvg(x,:) + posVecY(y)/specimen.height...
            *stress.xFirstMoment(x,:);
    end 
end

end

%{

% Get the iteration limits
[sizeY,sizeX,sizeT] = size(accel.x);

% Set the co-ord system central to the specimen
posVecX = pos.x;
posVecY = pos.y - specimen.height/2;
[posXGrid,posYGrid] = meshgrid(posVecX,posVecY);
posXGridF = padarray(posXGrid,[0,0,time.numFrames-1],'replicate','post'); 
posYGridF = padarray(posXGrid,[0,0,time.numFrames-1],'replicate','post'); 

% Calculate the weighted averages
ax_y = accel.x.*posYGridF;
ay_x = accel.y.*posXGridF;

% Pre-alloc for speed
avg_ax_y = zeros(sizeX,sizeT);
avg_ay_x = zeros(sizeX,sizeT);
avg_ay = zeros(sizeX,sizeT);

for t = 1:sizeT
    % Calculate the position weighted accleration fields for this frame
    ax_y(:,:,t) = accel.x(:,:,t).*posYGrid;
    ay_x(:,:,t) = accel.y(:,:,t).*posXGrid;
    
    for x = 1:sizeX
        % Calculate the surface average acceleration up to 'x'
        avg_ax_y(x,:) = mean(mean(ax_y(:,1:x,:)));
        avg_ay_x(x,:) = mean(mean(ay_x(:,1:x,:)));
        avg_ay(x,:) = mean(mean(accel.y(:,1:x,:)));
        
        % Calculate the first moment and the linear distribution of stress
        stress.xFirstMoment(x,:) = (12*material.rho*posVecX(x)/specimen.height)*...
            (avg_ax_y(x,:)-avg_ay_x(x,:)+posVecX(x)*avg_ay(x,:));
        
        for y = 1:sizeY
            stress.xLinearGauge(y,x,:) = stress.xAvg(x,:) + posVecY(y)/specimen.height...
                *stress.xFirstMoment(x,:);
        end 
    end
end
%}

