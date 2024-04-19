function data3D = func_reshapeAbaqusData(data,FEModel,sortInd,varCol)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 21/7/2017

    data3D = zeros(FEModel.numYNodes,FEModel.numXNodes,FEModel.numFrames);
    for f = 1:FEModel.numFrames
        temp = [];
        % Grab the frame and temporarily store it
        temp = data{f}(:,varCol);
        % Sort the data based on the node location
        temp = temp(sortInd,:);
        % Reshape the data into a frame
        temp = reshape(temp,[FEModel.numYNodes,FEModel.numXNodes]);
        % Save the data to a 3D matrix
        data3D(:,:,f) = temp;    
    end
end

