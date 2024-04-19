function [pos] = func_ExtrapolatePos(pos,ExPts)
%This function is written to pad the Pos Data Structure to match the size 
    %of the specimen

%  Author: Andrew Matejunas


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Raw Position
pos.Raw=pos;
frames=size(pos.Raw.xGridF,3);
%% Determine distance of extrapolation
dx=pos.xStep;
dy=pos.yStep;

%% Extrapolate x on left and right
xExLeftVec=flip(1:ExPts)*dx;
pos.x=pos.x;
xExLeft=pos.Raw.xGrid(:,1)-xExLeftVec;

xExRightVec=(1:ExPts)*dx;
xExRight=pos.Raw.xGrid(:,end)+xExRightVec;

pos.xGrid=[xExLeft,pos.Raw.xGrid,xExRight];

%% Extrapolate x on top and bottom
xExTop=repmat(pos.xGrid(1,:),ExPts,1);
xExBot=repmat(pos.xGrid(end,:),ExPts,1);
pos.xGrid=[xExTop;pos.xGrid;xExBot];

pos.x=squeeze(pos.xGrid(1,:));
pos.lengthX = pos.x(end)+pos.xStep/2;
pos.xGridF = padarray(pos.xGrid,[0,0,frames-1], ...
    'replicate','post');
pos.x0F = squeeze(padarray(pos.x,[0,0,frames-1], ...
    'replicate','post'));

%% Extrapolate Y on top and bottom
yExTopVec=flip(1:ExPts)'*dy;
yExTop=pos.Raw.xGrid(1,:)-yExTopVec;

yExBotVec=(1:ExPts)'*dy;
yExBot=pos.Raw.xGrid(end,:)+yExBotVec;
pos.yGrid=[yExTop;pos.Raw.yGrid;yExBot];
%% Extrapolate Y on left and right
yExLeft=repmat(pos.yGrid(:,1),1,ExPts);
yExRight=repmat(pos.yGrid(:,end),1,ExPts);
pos.yGrid=[yExLeft,pos.yGrid,yExRight];

pos.y=squeeze(pos.yGrid(:,1));
pos.yGridF = padarray(pos.yGrid,[0,0,frames-1], ...
    'replicate','post');

end