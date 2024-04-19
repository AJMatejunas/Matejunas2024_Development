function strengthFrameRange = func_manualInputFractureFrameRange(time)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 20/11/2017

    % Ask the user for the frame to obtain the strength
    while true
        frameData = inputdlg({'Strength Identification Frame (Min):','Strength Identification Frame (Max):'},...
            'Input frame range to obtain strength data', 1, {'1',num2str(time.numFrames)});
        strengthFrameMin = str2double(frameData{1});
        strengthFrameMax = str2double(frameData{2});
        % Make sure the user is giving a sensible input
        if ((strengthFrameMin<=time.numFrames) && (strengthFrameMin>=1)) && ...
            ((strengthFrameMax<=time.numFrames) && (strengthFrameMax>=1))
                break;
        else
            msg = 'Selected frame(s) do not exist!';
            h_mdlg = msgbox(msg,'Warning!');
            uiwait(h_mdlg);
        end
    end
    strengthFrameRange = strengthFrameMin:strengthFrameMax;
end

