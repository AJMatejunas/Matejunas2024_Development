function fracture = func_manualInputFractureLoc(pos,time)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 7/6/2017
%
% Prompts user to manually input the location and time of fracture
   
    % Ask the user for the frame to select the crack location 
    VarsToInput = {'Frame to Select Crack Location:','X Location (px):',...
        'Y Location (px):'};
    inputData = inputdlg(VarsToInput,'Input Fracture Data', 1,... 
        {num2str(time.numFrames),'1','1'});

    % Assign input variables
    fracture.locFrame = str2double(inputData{1});
    fracture.locX = str2double(inputData{2});
    fracture.locY = str2double(inputData{3});
end

