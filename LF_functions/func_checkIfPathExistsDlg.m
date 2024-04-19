function [plotSeq,existFlag] = func_checkIfPathExistsDlg(inPath,questStr,dlgTitle)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 19/2/2017

if exist(inPath,'file') == 7
    existFlag = true;
    choice = questdlg(questStr,dlgTitle,'Yes','No','No');
    switch choice
        case 'Yes'
            plotSeq = true;
        case 'No'
            plotSeq = false;
    end
else
    % If the directory does not exist, create it and plot the image seq
    existFlag = false;
    plotSeq = true;
    mkdir(inPath);
end

end

