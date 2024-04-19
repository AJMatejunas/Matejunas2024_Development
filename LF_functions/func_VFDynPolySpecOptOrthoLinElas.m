function [identStiffVF,VFDiagnostics] = func_VFDynPolySpecOptOrthoLinElas(VFOpts,...
    pos,specimen,material,accel,strain)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 18/9/2017
%
% Uses special optimised virtual fields to find the stiffness components
% Q11,Q22,Q12,Q66 for a dynamically impacted linear elastic orthotropic 
% material. Assumes the co-ordinate system is centred on the bottom LH 
% corner of the rectangular specimen such that all co-ords are positive. 
% Specimen has dimensions L*w and it is assumed the specimen is impacted 
% on the edge x1 = L.

% Turn off warnings for ill conditioned matrices
warning('off', 'MATLAB:singularMatrix'); 
warning('off', 'MATLAB:nearlySingularMatrix'); 

%--------------------------------------------------------------------------
% Parameter Initialisation
% Polynomial degrees
m = VFOpts.polyOrderX1;  % Polynomial order for x1
n = VFOpts.polyOrderX2;  % Polynomial order for x2
% Number of degrees of freedom for each virtual disp component 
nDof = (m+1)*(n+1);

% Number of speciality constraints, 2 for iso, 4 for ortho
numSpecConst = 4;

% Get the size of measurement array
[nRows,nCols,nFrames] = size(strain.x);
nPoints = nRows*nCols;

% Specimen Geometry
%L = specimen.length;
%w = specimen.height;
L = max(max(pos.xGrid));
w = max(max(pos.yGrid));
rho = material.rho;

% Vectorise the pos fields to speed up the calculation
x1 = reshape(pos.xGrid,[],1);
x2 = reshape(pos.yGrid,[],1);

%--------------------------------------------------------------------------
% Construct the Virtual Fields - Expand with Polynomials
% NOTE: 
% - This section is the same for iso and ortho models. 
% - Does not change in time, monomials are only a function of position.
% - BCs can be included in the polynomial by adding a factor e.g.:
%       To zero border at x1 = L use the factor (L-x1)
%       To zero border at x2 = w use the factor (w-x2)
%       This will require a re-calc of the derivatives for virtual strains

% Initialise to zero
epsStar1Mon = zeros(nPoints,nDof*2);
epsStar2Mon = epsStar1Mon; 
epsStar6Mon = epsStar1Mon;
uStar1Mon = epsStar1Mon;
uStar2Mon = epsStar1Mon;

% Loop over the polynomial orders to calculate monomial values
% NOTE: to keep dof order correct loop over j then i
dof = 1;
for j = 0:n
    for i = 0:m
        % Calculate the monomial derivative wrt x1
        dx1 = -1*(x1/L).^i.*(x2/w).^j + (L-x1).*i/L.*(x1/L).^(i-1).*(x2/w).^j;
        % Calculate the monomial derivative wrt x2
        dx2 = (L-x1).*(j/w*(x1/L).^i).*(x2/w).^(j-1);
        
        % Virtual strain fields - monomial componens
        epsStar1Mon(:,dof) = dx1;
        epsStar1Mon(:,dof+nDof) = zeros(nPoints,1);
        epsStar2Mon(:,dof) = zeros(nPoints,1);
        epsStar2Mon(:,dof+nDof) = dx2;
        epsStar6Mon(:,dof) = dx2;
        epsStar6Mon(:,dof+nDof) = dx1;
 
        % Virtual disp fields - monomial components
        % u1 only depends on Aij 
        uStar1Mon(:,dof) = (L-x1).*(x1/L).^i.*(x2/w).^j;
        uStar1Mon(:,dof+nDof) = zeros(nPoints,1);
        %u2 only depends on Bij
        uStar2Mon(:,dof) = zeros(nPoints,1);
        uStar2Mon(:,dof+nDof) = (L-x1).*(x1/L).^i.*(x2/w).^j;
        
        % Go to the next degree of freedom
        dof = dof+1;
    end
end

%--------------------------------------------------------------------------
% Construct the components of the H matrix - squares of the virtual strain
% comonents, NOTE: this does not change from iso to ortho
H11 = epsStar1Mon'*epsStar1Mon;
H22 = epsStar2Mon'*epsStar2Mon;
H12 = epsStar1Mon'*epsStar2Mon;
H66 = epsStar6Mon'*epsStar6Mon;

%----------------------------------------------------------------------
% Create Z vectors for solving the optimisation problem
% NOTE: only two conditions for an isotropic material
Aconst = [];
% Initialise as zeros
Z11=zeros(1,nDof*2+size(Aconst,1)+numSpecConst);
Z22=Z11; Z12=Z11; Z66=Z11;
% Speciality condition for Q11
Z11(end-numSpecConst+1:end) = [1,0,0,0];
% Speciality condition for Q22
Z22(end-numSpecConst+1:end) = [0,1,0,0];
% Speciality condition for Q12
Z12(end-numSpecConst+1:end) = [0,0,1,0];
% Speciality condition for Q66
Z66(end-numSpecConst+1:end) = [0,0,0,1];

%--------------------------------------------------------------------------
% Loop over each frame and solve the optimisation problem
startFrame = VFOpts.startFrame;

% Initialise stiffness vectors to zero 
Q11t = zeros(1,nFrames);
Q22t = Q11t; Q12t = Q11t; Q66t = Q11t; 
for ff = startFrame:nFrames % when testing only run for the first half
    % Vectorise the current frame of the strain and acceleration
    eps1 = reshape(squeeze(strain.x(:,:,ff)),[],1);
    eps2 = reshape(squeeze(strain.y(:,:,ff)),[],1);
    eps6 = reshape(squeeze(strain.s(:,:,ff)),[],1);
    accel1 = reshape(squeeze(accel.x(:,:,ff)),[],1);
    accel2 = reshape(squeeze(accel.y(:,:,ff)),[],1);
    
    %----------------------------------------------------------------------
    % Intialise the B matrices for speciality contraints - changes in time
    SC11 = zeros(1,nDof*2);
    SC22 = SC11; SC12 = SC11; SC66 = SC11;

    % Construct the vectors containing the speciality conditions 'SC'
    for kk = 1:nDof*2
       SC11(kk) = L*w*mean(epsStar1Mon(:,kk).*eps1);
       SC22(kk) = L*w*mean(epsStar2Mon(:,kk).*eps2);
       SC12(kk) = L*w*mean(epsStar2Mon(:,kk).*eps1+epsStar1Mon(:,kk).*eps2);
       SC66(kk) = L*w*mean(epsStar6Mon(:,kk).*eps6);
    end
    
    %----------------------------------------------------------------------
    % Assemble the constraint matrix - BCs and speciality
    A = [Aconst;SC11;SC22;SC12;SC66];
    B = zeros(size(A,1)); % Block of zeros for the optimisation matrix
    
    %//////////////////////////////////////////////////////////////////////
    % Solve the System Iteratively
    % NOTE: must be initialised to different vals for convergence
    Q = [10,1,0.5,0.5]*10^9;  % Initialise Q to initialise the first H matrix
    
    nIter = 20;         % max iters before breaking the loop
    deltaLim = 0.001;
    delta = 10;
    ii = 1;
    Qprev = Q;
    
    % Convergence loop
    while ii<nIter && delta>deltaLim
        % Assemble the Hessian matrix - Isotropic Case
        H = 2*(L*w/nPoints)^2*...
            ((Q(1)^2+Q(3)^2)*H11+...
            (Q(2)^2+Q(3)^2)*H22+...
            2*(Q(1)+Q(2))*Q(3)*H12+...
            Q(4)^2*H66);
       
        % Coefficient for normalising the matrices prior to inversion
        normCoeff = max(max(A))/max(max(H));

        % Assemble the M matrix to find the coeffs for the virtual fields using
        % lagrangian optimisation
        OptM = [H*normCoeff,A';A,B];

        % Solve for the poly coeffs Aij and Bij for each Qij virtual field
        Y11 = OptM\Z11';
        Y22 = OptM\Z22';
        Y12 = OptM\Z12';
        Y66 = OptM\Z66';

        % Remove the Lagrange multipliers
        Y11(nDof*2+1:end) = [];
        Y22(nDof*2+1:end) = [];
        Y12(nDof*2+1:end) = [];
        Y66(nDof*2+1:end) = [];
        
        % Calculate the virtual displacement field values at each point
        % using the identified poly coeffs and the monomials
        uStar1Eval11 = uStar1Mon*Y11;
        uStar2Eval11 = uStar2Mon*Y11;
        uStar1Eval22 = uStar1Mon*Y22;
        uStar2Eval22 = uStar2Mon*Y22;
        uStar1Eval12 = uStar1Mon*Y12;
        uStar2Eval12 = uStar2Mon*Y12;
        uStar1Eval66 = uStar1Mon*Y66;
        uStar2Eval66 = uStar2Mon*Y66;

        % Calculate the Qij values using the virtual disp fields and the
        % measured acceleration fields
        Q(1) = -rho*L*w*(mean(uStar1Eval11.*accel1) + mean(uStar2Eval11.*accel2));
        Q(2) = -rho*L*w*(mean(uStar1Eval22.*accel1) + mean(uStar2Eval22.*accel2));
        Q(3) = -rho*L*w*(mean(uStar1Eval12.*accel1) + mean(uStar2Eval12.*accel2));
        Q(4) = -rho*L*w*(mean(uStar1Eval66.*accel1) + mean(uStar2Eval66.*accel2));

        % Calculate the difference in the stiffnesses for convergence checking
        delta = sum((Qprev-Q).^2./Q.^2);
        
        % Move to the next iteration
        Qprev = Q;
        ii = ii + 1;
    end
    
    %----------------------------------------------------------------------
    % Store the diagnostics for the identification   
    nIterLog(ff) = ii; % Number of iterations for convergence
    % Hessian matrix for the final stiffnesses
    H = 2*(L*w/nPoints)^2*...
        ((Q(1)^2+Q(3)^2)*H11+...
        (Q(2)^2+Q(3)^2)*H22+...
        2*(Q(1)+Q(2))*Q(3)*H12+...
        Q(4)^2*H66);
    % Calculate the noise sensitivity
    eta11t(ff) = sqrt(Y11'*H*Y11);
    eta22t(ff) = sqrt(Y22'*H*Y22);
    eta12t(ff) = sqrt(Y12'*H*Y12);
    eta66t(ff) = sqrt(Y66'*H*Y66);
    
    %----------------------------------------------------------------------
    % Store the stiffness found for this frame
    Q11t(ff) = Q(1);
    Q22t(ff) = Q(2);
    Q12t(ff) = Q(3);
    Q66t(ff) = Q(4);
    
    % TODO: create function that plots the virtual fields for debugging
end

% Initialise the structure for returning the VF diagnostics
VFDiagnostics.nIterConv = nIterLog;
VFDiagnostics.eta11VsT = eta11t;
VFDiagnostics.eta22VsT = eta22t;
VFDiagnostics.eta12VsT = eta12t;
VFDiagnostics.eta66VsT = eta66t;

% Initialise the structure for returning the identified stiffness values
identStiffVF.Q11VsT = Q11t;
identStiffVF.Q22VsT = Q22t;
identStiffVF.Q12VsT = Q12t;
identStiffVF.Q66VsT = Q66t;

% Turn on warnings for ill conditioned matrices
warning('on', 'MATLAB:singularMatrix'); 
warning('on', 'MATLAB:nearlySingularMatrix');

end

