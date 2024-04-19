function strengthFrame = func_manualInputFractureFrame(pos,time)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 7/6/2017

    % Ask the user for the frame to obtain the strength
    while true
        frameData = inputdlg({'Input frame to obtain strength data:'}, ...
                 'Strength Identification Frame', 1, {num2str(time.numFrames)} );
        strengthFrame = str2double(frameData{1});
        % Make sure the user is giving a sensible input
        if (strengthFrame<=time.numFrames) && (strengthFrame>=1)
            break;
        else
            msg = 'Selected frame does not exist!';
            h_mdlg = msgbox(msg,'Warning!');
            uiwait(h_mdlg);
        end
    end

end

