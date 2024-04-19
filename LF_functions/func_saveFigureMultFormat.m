function func_saveFigureMultFormat(figHandle,saveFile,plotProps,plotParams)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 26/4/2019
% Date Edited: 14/6/2019
%
% Saves a figure in multiple formats, can use export_fig from the matlab
% online repository to save nicer vector graphics.

    print(figHandle,saveFile,plotProps.imageSeqFormat,plotProps.imageSeqSaveRes);
    if plotParams.saveImageMatFig
        saveas(figHandle,saveFile,'fig');
    end
    if plotParams.saveImageVecFig
        if exist('export_fig','file') == 2
            export_fig(saveFile,'-eps','-pdf','-c[Inf,Inf,Inf,Inf]',figHandle)
        else
            saveas(figHandle,saveFile,'pdf');
            saveas(figHandle,saveFile,'epsc');
        end
    end
end

