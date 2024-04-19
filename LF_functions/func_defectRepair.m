function [mapx,mapy] = func_defectRepair(mapx,mapy,N,m)
% Peform defect dectection on the displacement.  It compares the displacment
% value to each of its neighbors and when there is too much different, it
% replaces it by the average value of it neighbors. 
%
%Inputs:
% - mapx, map of displacement in the x direction
%
% - mapy, map of displacmeent in the y direction
%
% - N, number of points from the edge to ignore, creates a border around
%      the tests area
%
% - m, the selected standard deviation level.  The standard deviation is
% sorted for each time, t. The program uses the mth standard deviation as
% the error cut-off
%
%
% Ouputs:
% - mapx, corrected map of displacement in the x-direction
%
% - mapy, corrected map of displacement in the y-direction
%
%
% v.1 , from Haibin Zhu with modifications by Cedric Devivier for speed
% optimization - 13 March 2016
%
% v.2, changes by F. Davis, to make work inside current processing program - 21 March 2016 

%% Dim vars

[ny,nx,t_step] = size(mapx);

%% Begin Loop on time

for kk = 1:t_step
    
    % initializes variables
    f_x = zeros(ny, nx); f_y = zeros(ny, nx);
    g_x = zeros(ny, nx); g_y = zeros(ny, nx);
    
    
    %% Perform outlier detection
    
    %Set f_x and f_y equal to the displacement in the x and y directions
    %respectively with a border of N points around the map edge removed
    f_x(1+N:end-N,1+N:end-N) = squeeze(abs(mapx(1+N:end-N,1+N:end-N,t_step)));
    f_y(1+N:end-N,1+N:end-N) = squeeze(abs(mapy(1+N:end-N,1+N:end-N,t_step)));
    
    %Perform the 2D convolution of the displacement with a square pulse of
    %width 2*N + 1
    g_x(1+N:end-N,1+N:end-N) = conv2(abs(mapx),ones(2*N+1,2*N+1)/(2*N+1)^2,'valid');
    g_y(1+N:end-N,1+N:end-N) = conv2(abs(mapy),ones(2*N+1,2*N+1)/(2*N+1)^2,'valid');
    
    h_x = f_x-g_x;
    h_y = f_y-g_y;
    
    %Scale the measurement by the standard deviation
    s_x = abs((h_x-mean2(h_x))/std2(h_x)); % scaling
    s_y = abs((h_y-mean2(h_y))/std2(h_y)); % scaling
    
    %Set the threshold
    sx = sort(reshape(s_x,[],1),'descend');
    [pos_x(:,1),pos_x(:,2)] = find(s_x>sx(m+1));
    sy = sort(reshape(s_y,[],1),'descend');
    [pos_y(:,1),pos_y(:,2)] = find(s_y>sy(m+1));
    
    %% Replace the outliers with the mean value of the neighboring points
    
    %Loop and replace points in the x-displacement map
    m=size(pos_x,1);
    for w=1:m
        i=pos_x(w,1);
        j=pos_x(w,2);
        mapx(i,j,t_step)=(mapx(i-1,j,t_step)+mapx(i+1,j,t_step)+mapx(i,j+1,t_step)+mapx(i,j-1,t_step))/4;
    end
    
    %Loop and replace points in the y-displacement map
    n=size(pos_y,1);
    for w=1:n
        i=pos_y(w,1);
        j=pos_y(w,2);
        mapy(i,j,t_step)=(mapy(i-1,j,t_step)+mapy(i+1,j,t_step)+mapy(i,j+1,t_step)+mapy(i,j-1,t_step))/4;
    end
end

end
