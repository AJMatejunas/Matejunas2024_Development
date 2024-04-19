function [posSlice1,posSlice2] = func_initPosSliceStructs(specimen,pos,time,material,...
    angledSlice1,angledSlice2,debug)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 29/1/2019
% Date Edited: 29/1/2019
%
% Creates position structures for each slice which holds the moving
% co-ordinate systems required for the generalised stress-strain curve
% method for orthotropic elasticity.
    
    % Flag for plotting figures to check if the co-ords are correct
    if nargin < 7
        debug = false;
    end
    
    % Centre the co-ords on the centre of the specimen height    
    yCent = pos.y-specimen.height/2;
    
    % Only create the co-ords if we can actually take a slice without
    % intersecting the impact edge
    if angledSlice1.xMaxInd > 0
        
        posSlice1.mat11_0 = yCent/sind(material.rotAngle);
        posSlice1.mat11_0F = squeeze(padarray(posSlice1.mat11_0...
            ,[0,0,time.numFrames-1],'replicate','post')); 
        posSlice1.mat11_0GridF = squeeze(padarray(posSlice1.mat11_0'...
            ,[0,angledSlice1.xMaxInd-1,time.numFrames-1],'replicate','post')); 
        
        % Step along each slice and create an angled moving co-ordinate
        % system that follows along the sample length
        for ss = 1:angledSlice1.xMaxInd
            xShift = pos.x(ss)+(specimen.height/2)/tand(material.rotAngle);
            yShift = specimen.height/2;
            tempXGrid = pos.xGrid - xShift;
            tempYGrid = pos.yGrid - yShift;
            [posSlice1.mat11F{ss},posSlice1.mat22F{ss}] = func_rotateVector2D(...
                tempXGrid,tempYGrid,material.rotAngle);
        end
        
        % Calculate the co-ords of the slice in the global reference frame
        for xx = 1:angledSlice1.xMaxInd
            [posSlice1.globX(:,xx),posSlice1.globY(:,xx)] =...
                func_calcAngSliceCoords(xx,pos,angledSlice1);
        end
        posSlice1.globXF = squeeze(padarray(posSlice1.globX...
            ,[0,0,time.numFrames-1],'replicate','post')); 
        posSlice1.globYF = squeeze(padarray(posSlice1.globY...
            ,[0,0,time.numFrames-1],'replicate','post')); 
        
        if debug
            % Figures showing co-ord system to verify it moves with the
            % slice
            hf = figure; imagesc(posSlice1.mat11_0GridF(:,:,100)); colorbar; axis image; set(gca,'YDir','normal');
            hf = figure; imagesc(posSlice1.mat11F{1}); colorbar; axis image; set(gca,'YDir','normal'); title('Slice 1, 11, First Cut');
            hf = figure; imagesc(posSlice1.mat22F{1}); colorbar; axis image; set(gca,'YDir','normal'); title('Slice 1, 22, First Cut');
            hf = figure; imagesc(posSlice1.mat11F{angledSlice1.xMaxInd}); colorbar; axis image; set(gca,'YDir','normal'); title('Slice 1, 11, End Cut');
            hf = figure; imagesc(posSlice1.mat22F{angledSlice1.xMaxInd}); colorbar; axis image; set(gca,'YDir','normal'); title('Slice 1, 22, End Cut');
            posSlice1
        end
    else
        posSlice1 = nan;
    end
    
    % Only create the co-ords if we can actually take a slice without
    % intersecting the impact edge
    if angledSlice2.xMaxInd > 0
        posSlice2.mat22_0 = yCent/sind(90-material.rotAngle);
        posSlice2.mat22_0F = squeeze(padarray(posSlice2.mat22_0...
            ,[0,0,time.numFrames-1],'replicate','post')); 
        posSlice2.mat22_0GridF = squeeze(padarray(posSlice2.mat22_0'...
            ,[0,angledSlice2.xMaxInd-1,time.numFrames-1],'replicate','post')); 
        
        % Step along each slice and create an angled moving co-ordinate
        % system that follows along the sample length
        for ss = 1:angledSlice2.xMaxInd
            xShift = pos.x(ss)+(specimen.height/2)/tand(90-material.rotAngle);
            yShift = specimen.height/2;
            tempXGrid = pos.xGrid - xShift;
            tempYGrid = pos.yGrid - yShift;
            [posSlice2.mat11F{ss},posSlice2.mat22F{ss}] = func_rotateVector2D(...
                tempXGrid,tempYGrid,material.rotAngle);
        end
        
        % Calculate the co-ords of the slice in the global reference frame
        for xx = 1:angledSlice2.xMaxInd
            [posSlice2.globX(:,xx),posSlice2.globY(:,xx)] = ...
                func_calcAngSliceCoords(xx,pos,angledSlice2);
        end
        posSlice2.globXF = squeeze(padarray(posSlice2.globX...
            ,[0,0,time.numFrames-1],'replicate','post')); 
        posSlice2.globYF = squeeze(padarray(posSlice2.globY...
            ,[0,0,time.numFrames-1],'replicate','post')); 
        
        if debug
            % Figures showing co-ord system to verify it moves with the
            % slice
            hf = figure; imagesc(posSlice2.mat22_0GridF(:,:,100)); colorbar; axis image; set(gca,'YDir','normal');
            hf = figure; imagesc(posSlice2.mat11F{1}); colorbar; axis image; set(gca,'YDir','normal'); title('Slice 2, 11, First Cut');
            hf = figure; imagesc(posSlice2.mat22F{1}); colorbar; axis image; set(gca,'YDir','normal'); title('Slice 2, 22, First Cut');
            hf = figure; imagesc(posSlice2.mat11F{angledSlice2.xMaxInd}); colorbar; axis image; set(gca,'YDir','normal'); title('Slice 2, 11, End Cut');
            hf = figure; imagesc(posSlice2.mat22F{angledSlice2.xMaxInd}); colorbar; axis image; set(gca,'YDir','normal'); title('Slice 2, 22, End Cut');
            posSlice2
        end
    else
        posSlice2 = nan;
    end
end

