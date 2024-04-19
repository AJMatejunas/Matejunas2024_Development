function angledSlice = func_initAngledSliceStruct(specimen,pos,angledSlice)
% Author: L. Fletcher
% PhotoDyn Research Group
% Date Created: 11/12/2018 
% Date Edited: 11/12/2018

angledSlice.length = abs(specimen.height/sind(angledSlice.angle));
angledSlice.lengthX = abs(angledSlice.length*cosd(angledSlice.angle));
angledSlice.xMin = 0; 
angledSlice.xMax = specimen.length - angledSlice.lengthX;
angledSlice.xMinInd = 1;
%angledSlice.returnComps = [0,1,1]; 
angledSlice.length = abs(specimen.height/sind(angledSlice.angle));
angledSlice.lengthX = abs(angledSlice.length*cosd(angledSlice.angle));
angledSlice.xMax = specimen.length - angledSlice.lengthX;
angledSlice.xMaxInd = sum(pos.x < angledSlice.xMax);
angledSlice.xStep = pos.xStep;

end

