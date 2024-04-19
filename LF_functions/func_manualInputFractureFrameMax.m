function strengthFrameMax = func_manualInputFractureFrameMax(time)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 20/11/2017

    % Ask the user for the frame to obtain the strength
    while true
        frameData = inputdlg({'Strength Identification Frame (Max):'},...
            'Input frame range to obtain strength data', 1, {num2str(time.numFrames)});
        strengthFrameMax = str2double(frameData{1});
        % Make sure the user is giving a sensible input
        if ((strengthFrameMax<=time.numFrames) && (strengthFrameMax>=1))
                break;
        else
            msg = 'Selected frame(s) do not exist!';
            h_mdlg = msgbox(msg,'Warning!');
            uiwait(h_mdlg);
        end
    end
end

