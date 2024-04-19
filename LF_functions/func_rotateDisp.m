function [grid,disp] = func_rotateDisp(grid,disp,userInputFlag)   
% Author: Jared Van-Blitterswyk
% PhotoDyn Group, University of Southampton
% Date: 29/3/2017
% Edited By: Lloyd Fletcher 
    
    if nargin < 3
        userInputFlag = true;
    end
    
    if userInputFlag
        % Ask user to input the rotation angle
        while true
            inputData = inputdlg({'Rotation angle of inputData (degrees):'}, ...
                     'Input rotation angle', 1, {'0'} );
            grid.rotAngle = str2double(inputData{1});
            % Make sure the user is giving a sensible input
            if (grid.rotAngle >= -90) && (grid.rotAngle <= 90)
                break;
            else
                msg = 'Rotation angle needs to be between -90 and 90 degrees';
                h_mdlg = msgbox(msg,'Warning!');
                uiwait(h_mdlg);
            end
        end
    end
    
    % Rotate displacements
    % theta in degrees
    theta = grid.rotAngle;
    if theta ~= 0
        for i = 1:size(disp.x,1)
            for j = 1:size(disp.x,2)
                for k = 1:size(disp.x,3)
                    dispRot = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]*[disp.x(i,j,k); disp.y(i,j,k)];
                    disp.x(i,j,k) = dispRot(1);
                    disp.y(i,j,k) = dispRot(2);
                end
            end
        end
    end
end