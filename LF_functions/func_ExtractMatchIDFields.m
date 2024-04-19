function [disp,pos,DIC,specimen,strain] = func_ExtractMatchIDFields(DIC, ...
    specimen,material,globalOpts,ConvUnit)
%Extracts displacement, picture coordinates, and specimen parameters from
%Digital image Correlation data output by Match ID. 

%Author: Andrew Matejunas

%Function Inputs:
    %SS- Subset Size in Pixels
    %ST- Step Size in Pixels
    %specimen- Structure containing specimen geometric information
    %ConvUnit- Distance units of MatchID outputs to convert back to m
        %(e.g. 1e-3 for mm)
    

 %Outputs
    %pos- structure containing specimen coordinates
    %disp- structure containing x and y displacement arrays
    %DIC- structure containing DIC parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Choose File Locations
[xFile,xPath]=uigetfile('*', ...
    'Choose First File Containing X data coordinates');
[yFile,yPath]=uigetfile('*', ...
    'Choose First File Containing y data coordinates');

[uFile,uPath]=uigetfile('*', ...
    'Choose First File Containing X displacements');
[vFile,vPath]=uigetfile('*', ...
    'Choose First File Containing y displacements');

switch globalOpts.ExtractDICStrain
    case true
        [xxFile,xxPath]=uigetfile('*', ...
            'Choose First File Containing XX strains');
        [yyFile,yyPath]=uigetfile('*', ...
            'Choose First File Containing yy strains');
        [xyFile,xyPath]=uigetfile('*', ...
            'Choose First File Containing XY strains');
    case false
        fprintf('No strain extraction. Strains will be calculated from displacement fields \n')
end

%% Extract Image Coordinates
fprintf('Extracting x coordinates \n')
pos.xArray=func_readMatchIDMatrix(xPath,xFile);
%Match ID does not ouput coordinate system for the reference frame by
%default so here copy the first deformed image coordinates
pos.xArray(:,:,1)=pos.xArray(:,:,2)*ConvUnit;
pos.xGrid=squeeze(pos.xArray(:,:,1));
%fix Nans
for k=1:size(pos.xGrid,1)
pos.xGrid(k,:)=mean(pos.xGrid,"omitnan");
end
pos.x=squeeze(pos.xGrid(1,:));
pos.xStep=pos.x(2)-pos.x(1);

fprintf('Extracting y coordinates \n')
pos.yArray=func_readMatchIDMatrix(yPath,yFile);
%Match ID does not ouput coordinate system for the reference frame by
%default so here copy the first deformed image coordinates
pos.yArray(:,:,1)=pos.yArray(:,:,2)*ConvUnit;
pos.yGrid=squeeze(pos.yArray(:,:,1));
%fix Nans
for k=1:size(pos.yGrid,2)
pos.yGrid(:,k)=mean(pos.yGrid,2,"omitnan");
end
pos.y=squeeze(pos.yGrid(:,1));
pos.yStep=pos.y(2)-pos.y(1);


%% Extract Displacements
fprintf('Extracting x displacements \n')
disp.x=func_readMatchIDMatrix(uPath,uFile)*ConvUnit;

fprintf('Extracting y displacements \n')
disp.y=func_readMatchIDMatrix(vPath,vFile)*ConvUnit;

%% Extract strains
switch globalOpts.ExtractDICStrain
    case true
        fprintf('Extracting xx Strains \n')
        strain.x=func_readMatchIDMatrix(xxPath,xxFile);
        fprintf('Extracting yy Strains \n')
        strain.y=func_readMatchIDMatrix(yyPath,yyFile);
        fprintf('Extracting Shear Strains \n')
        strain.s=func_readMatchIDMatrix(xyPath,xyFile);
    case false
        strain=0;
end


%% Generate DIC Data Structure
DIC.numPtsX=length(pos.x);
DIC.numPtsY=length(pos.y);
DIC.mPerST=pos.xStep;
DIC.mPerPx=pos.xStep/DIC.ST;
DIC.SSpx=DIC.SS;
DIC.SSm=(DIC.ST/pos.xStep)*DIC.SS;
DIC.STpx=DIC.ST;
DIC.STm=pos.xStep;

%% Update Specimen Geometry
fprintf('Updating Speciment Geometry \n')
[specimen,DIC]=func_updateSpecGeomDIC_v4(specimen,material,DIC,disp);

end