function [identStiffVF,VFDiagnostics] = func_VFDynPolySpecOptIsoLinElas(VFOpts,...
    pos,specimen,material,accel,strain)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 18/9/2017
%
% Uses special optimised virtual fields to find the stiffness components
% Q11 and Q12 for a dynamically impacted linear elastic isotropic material. 
% Assumes the co-ordinate system is centred on the bottom LH corner of the
% rectangular specimen such that all co-ords are positive. Specimen has
% dimensions L*w and it is assumed the specimen is impacted on the edge x1
% = L.

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
numSpecConst = 2;

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
Z12=Z11;
% Speciality condition for Q11
Z11(end-numSpecConst+1:end) = [1,0];
% Speciality condition for Q12
Z12(end-numSpecConst+1:end) = [0,1];

%--------------------------------------------------------------------------
% Loop over each frame and solve the optimisation problem
startFrame = VFOpts.startFrame;

% Initialise stiffness vectors to zero 
Q11t = zeros(1,nFrames);
Q12t = Q11t; Et = Q11t; Nut = Q11t; 
for ff = startFrame:nFrames % when testing only run for the first half
    % Vectorise the current frame of the strain and acceleration
    eps1 = reshape(squeeze(strain.x(:,:,ff)),[],1);
    eps2 = reshape(squeeze(strain.y(:,:,ff)),[],1);
    eps6 = reshape(squeeze(strain.s(:,:,ff)),[],1);
    accel1 = reshape(squeeze(accel.x(:,:,ff)),[],1);
    accel2 = reshape(squeeze(accel.y(:,:,ff)),[],1);
    
    %----------------------------------------------------------------------
    % Intialise the B matrices for speciality contraints - changes in time
    % NOTE: two constraints for isotropic materials
    SC11 = zeros(1,nDof*2);
    SC12 = zeros(1,nDof*2);

    % Construct the components of the B matrix
    for kk = 1:nDof*2
       SC11(kk) = L*w*mean(epsStar1Mon(:,kk).*eps1 + epsStar2Mon(:,kk).*eps2 + ...
           0.5*epsStar6Mon(:,kk).*eps6);
       SC12(kk) = L*w*mean(epsStar2Mon(:,kk).*eps1 + epsStar1Mon(:,kk).*eps2 - ...
           0.5*epsStar6Mon(:,kk).*eps6);
    end

    %----------------------------------------------------------------------
    % Assemble the constraint matrix - BCs and speciality
    A = [Aconst;SC11;SC12];
    B = zeros(size(A,1)); % Block of zeros for the optimisation matrix
    
    %//////////////////////////////////////////////////////////////////////
    % Solve the System Iteratively
    % NOTE: must be initialised to different vals for convergence
    Q = [2e9,1e9];  % Initialise Q to initialise the first H matrix
    
    nIter = 20;         % max iters before breaking the loop
    deltaLim = 0.001;
    delta = 10;
    ii = 1;
    Qprev = Q;
    
    % Convergence loop
    while ii<nIter && delta>deltaLim
        % Assemble the Hessian matrix - Isotropic Case
        H = (L*w/nPoints)^2*(...
           (Q(1)^2+Q(2)^2)*H11+...
           (Q(1)^2+Q(2)^2)*H22+...
           4*Q(1)*Q(2)*H12+...
           ((Q(1)-Q(2))/2)^2*H66);
       
        % Coefficient for normalising the matrices prior to inversion
        normCoeff = max(max(A))/max(max(H));

        % Assemble the M matrix to find the coeffs for the virtual fields using
        % lagrangian optimisation
        OptM = [H*normCoeff,A';A,B];

        % Solve for the poly coeffs Aij and Bij for each Qij virtual field
        Y11 = OptM\Z11';
        Y12 = OptM\Z12';

        % Remove the Lagrange multipliers
        Y11(nDof*2+1:end) = [];
        Y12(nDof*2+1:end) = [];
        
        % Calculate the virtual displacement field values at each point
        % using the identified poly coeffs and the monomials
        uStar1Eval1 = uStar1Mon*Y11;
        uStar2Eval1 = uStar2Mon*Y11;
        uStar1Eval2 = uStar1Mon*Y12;
        uStar2Eval2 = uStar2Mon*Y12;

        % Calculate the Qij values using the virtual disp fields and the
        % measured acceleration fields
        Q(1) = -rho*L*w*(mean(uStar1Eval1.*accel1) + mean(uStar2Eval1.*accel2));
        Q(2) = -rho*L*w*(mean(uStar1Eval2.*accel1) + mean(uStar2Eval2.*accel2));

        % Calculate the difference in the stiffnesses for convergence checking
        delta = sum((Qprev-Q).^2./Q.^2);
        
        % Move to the next iteration
        Qprev = Q;
        ii = ii + 1;
    end
    
    %----------------------------------------------------------------------
    % Store the diagnostics for the identification
    % Number of iterations for convergence
    nIterLog(ff) = ii;
    H = (L*w/nPoints)^2*(...
       (Q(1)^2+Q(2)^2)*H11+...
       (Q(1)^2+Q(2)^2)*H22+...
       4*Q(1)*Q(2)*H12+...
       ((Q(1)-Q(2))/2)^2*H66);
    % Sensitivity of each identified parameter
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

% Initialise the structure for returning the VF diagnostics
VFDiagnostics.nIterConv = nIterLog;
VFDiagnostics.etaxxVsT = eta11t;
VFDiagnostics.etaxyVsT = eta12t;

% Initialise the structure for returning the identified stiffness values
identStiffVF.QxxVsT = Q11t;
identStiffVF.QxyVsT = Q12t;
identStiffVF.ExxVsT = Et;
identStiffVF.NuxyVsT = Nut;

% Turn on warnings for ill conditioned matrices
warning('on', 'MATLAB:singularMatrix'); 
warning('on', 'MATLAB:nearlySingularMatrix'); 

end

