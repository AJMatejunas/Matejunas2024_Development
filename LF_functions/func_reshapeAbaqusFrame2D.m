function data2D = func_reshapeAbaqusFrame2D(data,numNodes,sortInd,varCol)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 9/8/2017

    % Grab the frame and temporarily store it
    temp = data{1}(:,varCol);
    % Sort the data based on the node location
    temp = temp(sortInd,:);
    % Reshape the data into a frame
    temp = reshape(temp,[numNodes.y,numNodes.x]);
    % Save the data to a 3D matrix
    data2D = temp;    
end


