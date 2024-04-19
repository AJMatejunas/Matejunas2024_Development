function [data,nodeNums] = func_loadAbaqusReport(reportFile,fileParams)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 21/7/2017
%
% Reads abaqus report file based on number of nodes and ignores repeated
% headers in the file based on the parameters passed in the fileParams
% struct

% Based on the number of variables we are loading create the search string
scanStr = '%f ';
scanStr = repmat(scanStr,[1,fileParams.numVars]);

f = 1;
fHandle = fopen(reportFile);
while true
    %fprintf('Reading frame: %i\n',f)

    if f == 1
        tempData = textscan(fHandle,scanStr,fileParams.numDataPts,...
            'headerLines',(fileParams.headerRows-1));
    else
        tempData = textscan(fHandle,scanStr,fileParams.numDataPts,...
            'headerLines',(fileParams.headerRows+fileParams.padRows));
    end
    
    if isempty(tempData{1})
        break;
    else
        nodeNums = tempData{1};
        
        dataMatrix = [];
        for i = 2:length(tempData)
            dataMatrix = [dataMatrix,tempData{i}];
        end
        data{f} = dataMatrix;
        f = f+1;
    end
end
fclose(fHandle);

end

