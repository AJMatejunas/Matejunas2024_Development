function [identStiffVFOpt,VFOptDiag] = func_VFDynPWSpecOptIsoLinElas(VFOpts,...
    pos,specimen,material,accel,strain)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date Created: 19/9/2017
% Date Edited: 5/8/2019
%
% Uses piece wise, bi-linear FE special optimised virtual fields to find 
% the stiffness components Q11,Q12 for a dynamically impacted linear 
% elastic isotropic material. Specimen has dimensions L*w and it is assumed 
% the specimen is impacted on the edge x1 = L.

% Remove the extrapolated data on the impact edge if needed
if VFOpts.cutImpactEdge
    [pos,accel,strain] = func_removeFieldImpactEdge(VFOpts,pos,accel,strain);
end

% Check if we are blocking just the impact edge or both
if ~isfield(VFOpts,'numBCs')
    VFOpts.numBCs = 1;
end

% Turn off warnings for ill conditioned matrices
warning('off', 'MATLAB:singularMatrix'); 
warning('off', 'MATLAB:nearlySingularMatrix');

%--------------------------------------------------------------------------
% Parameter Initialisation
% Number of elements in the virtual mesh
m = VFOpts.nElemX;  % number of elements along x1
n = VFOpts.nElemY;  % number of elements along x2

% Geometry and material properties
rho = material.rho;
L = pos.x(end)+pos.xStep/2;
w = pos.y(end)+pos.yStep/2;
% NOTE: L and w have to use pos if the impact edge is cropped!

% Number of speciality constraints, 2 for isotropic
numSpecConst = 2;

% Mesh geometry parameters
nNodes = (m+1)*(n+1);
nElems = m*n;
[nRow,nCol,nFrames] = size(strain.x); 
nPoints = nCol*nRow;
lElem = L/m;
wElem = w/n;

% Vectorise the fields to speed up the calculation
% NOTE: in the FE formulation co-ord system is critical, co-ords should be
% all be positive
x1 = reshape(pos.xGrid,nPoints,1);
x2 = reshape(pos.yGrid,nPoints,1);

%--------------------------------------------------------------------------
% Construct the Virtual Fields
% Piecewise functions - Bi-linear 4-node elements
% Element col number for element along x1
elemColNum = floor(x1*m/L)+1;
% Element row number for element along x2
elemRowNum = floor(x2*n/w)+1;

% Parametric co-ords in element co-ord system
x1Para = 2*x1/lElem-elemColNum*2+1;
x2Para = 2*x2/wElem-elemRowNum*2+1;

% Virtual Displacement Calculation
u0 = zeros(nPoints,1);
u1StarElem = 0.25*[(1-x1Para).*(1-x2Para) , u0,...
                   (1+x1Para).*(1-x2Para) , u0,...
                   (1+x1Para).*(1+x2Para) , u0,...
                   (1-x1Para).*(1+x2Para) , u0];
u2StarElem = 0.25*[u0, (1-x1Para).*(1-x2Para),...
                   u0, (1+x1Para).*(1-x2Para),...
                   u0, (1+x1Para).*(1+x2Para),...
                   u0, (1-x1Para).*(1+x2Para)];

% Virtual Strain Calculation
v0 = zeros(nPoints,1);
eps1StarElem = 1/2/lElem*[-(1-x2Para), v0,...
                           (1-x2Para), v0,...
                           (1+x2Para), v0,...
                          -(1+x2Para), v0];
eps2StarElem = 1/2/wElem*[v0, -(1-x1Para),...
                          v0, -(1+x1Para),...
                          v0,  (1+x1Para),...
                          v0,  (1-x1Para)];
eps6StarElem = 1/2*[-(1-x1Para)/wElem, -(1-x2Para)/lElem,...
                    -(1+x1Para)/wElem,  (1-x2Para)/lElem,...
                     (1+x1Para)/wElem,  (1+x2Para)/lElem,...
                     (1-x1Para)/wElem, -(1+x2Para)/lElem];

%--------------------------------------------------------------------------
% Construct Speciality Conditions and Components of the Hessian Matrix
% SCij are the speciality conditions
% Hij and the components of the Hessian

% Intitalise the SCij and Hij vars to zero
SC11 = zeros(1,2*nNodes); SC12 = SC11;
H11 = zeros(2*nNodes,2*nNodes);
H22 = H11; H12 = H11; H66 = H11;

% Definition of Node Connectivity
node1 = (elemColNum-1)*(n+1)+elemRowNum;     % First node of each elem (bottom LH)
node2 = elemColNum*(n+1)+elemRowNum;         % Second node of each elem (bottom RH)
node3 = elemColNum*(n+1)+elemRowNum+1;       % Third node of each elem (top RH)
node4 = (elemColNum-1)*(n+1)+elemRowNum+1;   % Fourth node of each elem (top LH)
% Matrix giving dofs affected by each data point
assemMat = [2*node1-1, 2*node1, 2*node2-1, 2*node2,...
            2*node3-1, 2*node3, 2*node4-1, 2*node4];

% Build the nodal components of the Hessian matrix, Hij 
for pp = 1:nPoints
    assemPt = assemMat(pp,:);   
    H11(assemPt,assemPt) = H11(assemPt,assemPt) + ...
        eps1StarElem(pp,:)'* eps1StarElem(pp,:);
    H22(assemPt,assemPt) = H22(assemPt,assemPt) + ...
        eps2StarElem(pp,:)'* eps2StarElem(pp,:);
    H12(assemPt,assemPt) = H12(assemPt,assemPt) + ...
        eps1StarElem(pp,:)'* eps2StarElem(pp,:);
    H66(assemPt,assemPt) = H66(assemPt,assemPt) + ...
        eps6StarElem(pp,:)'* eps6StarElem(pp,:); 
end

%--------------------------------------------------------------------------
% Define the Virtual Boundary Conditions and Constraint Matrix
% Each row of the a matrx is a linear constraint equation
Aconst = zeros(VFOpts.numBCs*(n+1),2*nNodes);

% u1*=0 on the right boundary (impacted edge of the sample)
for bb = 1:(n+1)
    Aconst(bb,2*nNodes-2*(n+1)+2*bb-1)=1;
end

if VFOpts.numBCs == 2
    % u1*=0 on the left boundary (free edge of the sample)
    for bb = 1:(n+1)
        % Start after the first n+1 equations
        Aconst(bb+(n+1),2*bb-1)=1;
    end
end

% Block of zeros to fill out M matrix
B = zeros(size(Aconst,1)+numSpecConst); 

%--------------------------------------------------------------------------
% Construct the Z vectors for solving the system
Z11 = zeros(1,2*nNodes+size(Aconst,1)+numSpecConst); Z12=Z11; 
Z11(end-numSpecConst+1:end) = [1,0];      % Speciality condition for Q11
Z12(end-numSpecConst+1:end) = [0,1];      % Speciality condition for Q12

%--------------------------------------------------------------------------
% Initialise the loop storage vars to zero
Q11t = zeros(1,nFrames);
Q12t = zeros(1,nFrames);
Et = zeros(1,nFrames); 
Nut = zeros(1,nFrames);
nIterLog = zeros(1,nFrames); 
mMatrixCond = zeros(1,nFrames); 
eta11t = zeros(1,nFrames);
eta12t = zeros(1,nFrames);

%--------------------------------------------------------------------------
% Loop over each frame and identify the stiffness
startFrame = VFOpts.startFrame;
for ff = startFrame:VFOpts.endFrame
    %----------------------------------------------------------------------
    % Reshape the strain and accel fields for this frame
    eps1 = reshape(strain.x(:,:,ff),nPoints,1);
    eps2 = reshape(strain.y(:,:,ff),nPoints,1);
    eps6 = reshape(strain.s(:,:,ff),nPoints,1);
    accel1 = reshape(accel.x(:,:,ff),nPoints,1);
    accel2 = reshape(accel.y(:,:,ff),nPoints,1);
    
    %----------------------------------------------------------------------
    % Assemble the speciality conditions and constraint matrix 
    
    % Set initial SCij vector to zero
    SC11 = zeros(1,2*nNodes); SC12 = SC11;
    % Calculate SCij summing over all points
    for pp = 1:nPoints
        assemPt = assemMat(pp,:);
        SC11(assemPt) = SC11(assemPt)+(L*w/nPoints)*...
            (eps1(pp)*eps1StarElem(pp,:)+eps2(pp)*eps2StarElem(pp,:)...
            +0.5*eps6(pp)*eps6StarElem(pp,:));
        SC12(assemPt) = SC12(assemPt)+(L*w/nPoints)*...
            (eps1(pp)*eps2StarElem(pp,:)+eps2(pp)*eps1StarElem(pp,:)...
            -0.5*eps6(pp)*eps6StarElem(pp,:));
    end
    
    % Combine all constraints into the A matrix, BCs and speciality
    A = [Aconst; SC11; SC12];

    %----------------------------------------------------------------------
    % Solve the System Iteratively
    Q = [2,1]*10^9;     % Initialise Q to initialise the first H matrix
    maxIter = 10;       % Max number of iters before breaking loop
    deltaLim = 0.001;   % Tolerance for convergence
    delta = 10;         % Set initial difference value to higher than tol
    Qprev = Q;          % Intialise prev Q matrix for conv check
    
    % Start Convergence loop
    ii = 1;
    while ii<maxIter && delta>deltaLim
        % Assemble the Hessian matrix - Isotropic Case
        H = (L*w/nPoints)^2*(...
           (Q(1)^2+Q(2)^2)*H11+...
           (Q(1)^2+Q(2)^2)*H22+...
           4*Q(1)*Q(2)*H12+...
           ((Q(1)-Q(2))/2)^2*H66);
        
        % Coefficient for normalising the matrices prior to inversion
        normCoeff = max(max(A))/max(max(H));

        % Assemble the M matrix to find the nodal displacements for the
        % piecewise virtual fields using lagrange optimisation
        OptM = [H*normCoeff,A';A,B];
        
        % Solve for the virtual nodal displacements for each Qij virtual field
        Y11 = OptM\Z11';
        Y12 = OptM\Z12';

        % Remove the Lagrange multipliers
        Y11(2*nNodes+1:end) = [];
        Y12(2*nNodes+1:end) = [];
        
        % Pre-alloc for speed and clear the vars to zero
        uStar1Eval1 = zeros(1,nPoints);
        uStar2Eval1 = uStar1Eval1;
        uStar1Eval2 = uStar1Eval1;
        uStar2Eval2 = uStar1Eval1;
        % Evaluate the optimised virtual displacement fields at all
        % measurement points
        for pp = 1:nPoints
            u1vv = zeros(1,2*nNodes);
            u2vv = zeros(1,2*nNodes);
            assemPt = assemMat(pp,:);
            u1vv(assemPt) = u1StarElem(pp,:);
            u2vv(assemPt) = u2StarElem(pp,:);
            % Virtual disp fields for Q11
            uStar1Eval1(pp) = u1vv*Y11;
            uStar2Eval1(pp) = u2vv*Y11;
            % Virtual disp fields for Q12
            uStar1Eval2(pp) = u1vv*Y12;
            uStar2Eval2(pp) = u2vv*Y12;
        end
        
        % Calculate the Qij values using the virtual disp fields and the
        % measured acceleration fields
        Q(1) = -rho*L*w*(mean(uStar1Eval1'.*accel1) + mean(uStar2Eval1'.*accel2));
        Q(2) = -rho*L*w*(mean(uStar1Eval2'.*accel1) + mean(uStar2Eval2'.*accel2));

        % Calculate the difference in the stiffnesses for convergence checking
        delta = sum((Qprev-Q).^2./Q.^2);
        
        % Move to the next iteration
        Qprev = Q;
        ii = ii + 1;
    end
    %----------------------------------------------------------------------
    % Store diagnostics for the identification
    nIterLog(ff) = ii;
    mMatrixCond(ff) = cond(OptM);
    % Hessian matrix for the final stiffnesses
    H = (L*w/nPoints)^2*(...
       (Q(1)^2+Q(2)^2)*H11+...
       (Q(1)^2+Q(2)^2)*H22+...
       4*Q(1)*Q(2)*H12+...
       ((Q(1)-Q(2))/2)^2*H66);
    % Calculate the noise sensitivity
    eta11t(ff) = sqrt(Y11'*H*Y11);
    eta12t(ff) = sqrt(Y12'*H*Y12);
    
    %----------------------------------------------------------------------
    % Store the stiffnesses found for this frame
    Q11t(ff) = Q(1);
    Q12t(ff) = Q(2); 
    % Convert these to E and nu
    Nut(ff) = Q(2)/Q(1);
    Et(ff) = Q(1)*(1-Nut(ff)^2);
end

% Push the diagnostics into a data struct to return
VFOptDiag.nIters = nIterLog;
VFOptDiag.mMatrixCond = mMatrixCond; 
VFOptDiag.eta11VsT = eta11t;
VFOptDiag.eta12VsT = eta12t;

% Push the identified stiffnesses into a struct to return
identStiffVFOpt.QxxVsT = Q11t;
identStiffVFOpt.QxyVsT = Q12t;
identStiffVFOpt.ExxVsT = Et;
identStiffVFOpt.NuxyVsT = Nut;

% Calculate medians over the identification time to reject outliers
% Note: this is a first guess at the identified value
identStiffVFOpt.ExxAvgOverT = nanmedian(identStiffVFOpt.ExxVsT);
identStiffVFOpt.NuxyAvgOverT = nanmedian(identStiffVFOpt.NuxyVsT);
identStiffVFOpt.QxxAvgOverT = nanmedian(identStiffVFOpt.QxxVsT);
identStiffVFOpt.QxyAvgOverT = nanmedian(identStiffVFOpt.QxyVsT);

% Turn warnings back on for ill conditioned matrices
warning('on', 'MATLAB:singularMatrix'); 
warning('on', 'MATLAB:nearlySingularMatrix');

end

