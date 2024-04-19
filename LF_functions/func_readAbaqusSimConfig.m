function FEModel = func_readAbaqusSimConfig(configFile,fileParams,matModel,loading)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 21/7/2017
% Date Edited: 22/1/2019
%
% Reads user generated abaqus simulation configuration (.txt) file and 
% returns the data in the FEModel struct.


    fHandle = fopen(configFile);
    
    timeData = textscan(fHandle,'%f %f %f','headerLines',fileParams.headerRows);
    frameData = textscan(fHandle,'%f %f %f %f','headerLines',fileParams.headerRows-1);
    dampingData = textscan(fHandle,'%f %f %f','headerLines',fileParams.headerRows-1);
    FEModel.time.total = timeData{1};
    FEModel.time.solvStep = timeData{2};
    FEModel.time.outputStep = timeData{3};
    FEModel.frames.total = frameData{1};
    FEModel.frames.test = frameData{2};
    FEModel.frames.waveGuide = frameData{3};
    FEModel.frames.blank = frameData{4};
    FEModel.damping.beta = dampingData{1};
    FEModel.damping.bVisL = dampingData{2};
    FEModel.damping.bVisQ = dampingData{3};
    
    if strcmp(matModel,'ortho')
        specGeom = textscan(fHandle,'%f %f %f %f','headerLines',fileParams.headerRows);
        specMat = textscan(fHandle,'%f %f %f %f %f %f','headerLines',fileParams.headerRows-1);
        FEModel.spec.elemSize = specGeom{1};
        FEModel.spec.length = specGeom{2};
        FEModel.spec.height = specGeom{3};
        FEModel.spec.thickness = specGeom{4};
        FEModel.spec.rho = specMat{1};
        FEModel.spec.E11 = specMat{2};
        FEModel.spec.nu12 = specMat{3};
        FEModel.spec.E22 = specMat{4};
        FEModel.spec.G12 = specMat{5};
        FEModel.spec.rotAngle = specMat{6};
    else
        specGeom = textscan(fHandle,'%f %f %f %f','headerLines',fileParams.headerRows);
        specMat = textscan(fHandle,'%f %f %f','headerLines',fileParams.headerRows-1);
        FEModel.spec.elemSize = specGeom{1};
        FEModel.spec.length = specGeom{2};
        FEModel.spec.height = specGeom{3};
        FEModel.spec.thickness = specGeom{4};
        FEModel.spec.rho = specMat{1};
        FEModel.spec.E = specMat{2};
        FEModel.spec.nu = specMat{3};
    end

    if strcmp(loading,'projectile')
        projGeom = textscan(fHandle,'%f %f %f %f','headerLines',fileParams.headerRows);
        projMat = textscan(fHandle,'%f %f %f %f','headerLines',fileParams.headerRows-1);
        FEModel.proj.elemSize = projGeom{1};
        FEModel.proj.length = projGeom{2};
        FEModel.proj.height = projGeom{3};
        FEModel.proj.thickness = projGeom{4};
        FEModel.proj.vel = projMat{1};
        FEModel.proj.rho = projMat{2};
        FEModel.proj.E = projMat{3};
        FEModel.proj.nu = projMat{4};

        guideGeom = textscan(fHandle,'%f %f %f %f','headerLines',fileParams.headerRows);
        guideMat = textscan(fHandle,'%f %f %f','headerLines',fileParams.headerRows-1);
        FEModel.guide.elemSize = guideGeom{1};
        FEModel.guide.length = guideGeom{2};
        FEModel.guide.height = guideGeom{3};
        FEModel.guide.thickness = guideGeom{4};
        FEModel.guide.rho = guideMat{1};
        FEModel.guide.E = guideMat{2};
        FEModel.guide.nu = guideMat{3};
    end

    fclose(fHandle);

end

