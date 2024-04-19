function [pos,disp,strain] = func_loadDispStrainData(path,file)
% Loads grid GUI style matlab file into data structs for processing
load([path,file]);

% Push the Image Def grid data into structs
pos.x = X(1,:)*10^-3;
pos.xStep = pos.x(2) - pos.x(1); 
pos.y = Y(:,1)'*10^-3;
pos.yStep = pos.y(2) - pos.y(1); 
disp.x = disp_x*10^-3;
disp.y = disp_y*10^-3;

strain.x = strain_x*10^-3;
strain.y = strain_y*10^-3;
strain.s = strain_s*10^-3;

end

