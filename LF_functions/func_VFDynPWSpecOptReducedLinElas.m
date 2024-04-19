function [identStiffVFOpt,VFOptDiag,VFs] = func_VFDynPWSpecOptReducedLinElas(VFOpts, pos, strain, accel, specimen, material)
% Author: Jared Van-Blitterswyk
% PhotoDyn Group, University of Southampton
% Date Edited: 21/9/2018
%
% Piecewise optimised virtual fields method for identifying a single
% stiffness parameter (e.g. transverse or intralaminar composite).

% Remove the extrapolated data on the impact edge if needed
if VFOpts.cutImpactEdge
    [pos,accel,strain] = func_removeFieldImpactEdge(VFOpts,pos,accel,strain);
end

% Turn off warnings for ill conditioned matrices
warning('off', 'MATLAB:singularMatrix'); 
warning('off', 'MATLAB:nearlySingularMatrix');

% choose the number of elements in the x and y direction
m = VFOpts.nElemX; % number of elements in the X direction
n = VFOpts.nElemY; % number of elements in the Y direction
L = pos.x(end)+pos.xStep/2;
w = pos.y(end)+pos.yStep/2;
[rows, cols, frames] = size(strain.x);

n_nodes = (m+1)*(n+1); % number of nodes
n_points = rows*cols; % number of measurement points

L_el = L/m; % length of elements
w_el = w/n; % width of elements

%reshape coordinate matrices into vectors
X1 = reshape(pos.xGrid, n_points,1);
X2 = reshape(flipud(pos.yGrid), n_points,1);

%% Build Up Virtual Fields
% assign an element to each measurement point
iii = floor(X1*m/L)+1; % elements in the x1 direction
jjj = floor(X2*n/w)+1; % elements in the x2 direction

%define parametric coordinates
xsi1 = 2*X1/L_el - 2*iii +1; % along x1 direction
xsi2 = 2*X2/w_el - 2*jjj +1; % along x2 direction

% calculate virtual displacements - 1 dof per node
u1elem = 0.25*[(1-xsi1).*(1-xsi2) (1+xsi1).*(1-xsi2) (1+xsi1).*(1+xsi2) (1-xsi1).*(1+xsi2)]; % first row of N matrix (condensed version from pg. 51)
    
% calculate virtual strains
Eps1elem = [-(1-xsi2) (1-xsi2) (1+xsi2) -(1+xsi2)]*1/2/L_el;

%% Consruct the matrices in the form of the final optimization matrix 
% Bij will be used for the speciality condition
% Hij will be used to compute the Hessian matrix
B11 = zeros(1,n_nodes);
H11 = zeros(n_nodes, n_nodes);

% Define the nodes
n1 = (iii-1)*(n+1) + jjj; % first node
n2 = iii*(n+1) + jjj; % second node
n3 = iii*(n+1) + jjj +1; % third node
n4 = (iii-1)*(n+1) + jjj +1; % fourth node

% matrix containing the degrees of freedom affected by each data point
assemble = [n1 n2 n3 n4];

%% Define Virtual Boundary Conditions
Aconst = zeros(m*n+n+1, n_nodes); %there are (n+1) boundary conditions and 1 degrees of freedom per node
% for i = 1:(n+1) % u1(x1 = 0) = 0
%     Aconst(i,i) = 1; % set to one so that Aconst*Y yields equation describing constraints
% end
for i = 1:(n+1) % u1(x1 = L) = 0 (n+1) conditions
    Aconst(i, n_nodes-(n+1)+i) = 1;
end

% constrain horizontal virtual displacements to be the same within each
% vertical plane
for j = 1:m % logic written up on pg. 134 of logbook
    for i = 1:n
        Aconst(i+n+1+(j-1)*n, n_nodes-(m+2-j)*(n+1)+i) = 1;
        Aconst(i+n+1+(j-1)*n, n_nodes-(m+2-j)*(n+1)+i+1) = -1;
    end
end

%% Construct Z vector of zeros except for conditions of speciality
% Za is the vector used to find the virtual field to identify Q11
Za = zeros(1, n_nodes + size(Aconst,1)+1); % 32 columns plus 16 for describing the boundary conditions
Za(n_nodes + size(Aconst,1)+1 : n_nodes + size(Aconst,1)+1) = 1; % special VF condition for Q11

%% Compute stiffness for each frame using VFs with minimal noise sensitivity

% Allocate memory for vector holding stiffness values as a function of time
nIterLog = zeros(1,frames);
mMatrixCond = zeros(1,frames);
Q11t = zeros(1,frames);
u1va = zeros(1,n_points);

for j = VFOpts.startFrame:VFOpts.endFrame
    % reshape matrices for each frame9:end-9,1:end-9,:
    Eps1 = reshape(squeeze(strain.x(:,:,j)), n_points,1);
    A1 = reshape(squeeze(accel.x(:,:,j)), n_points,1);

    for k = 1:n_points
        assemble1 = assemble(k,:);
        B11(assemble1) = B11(assemble1) + (Eps1(k)*Eps1elem(k,:))*(L*w/n_points); % multiply actual strains by virtual strains

        % assemble Hessian matrix (minimization of sensitivity to noise)
        H11(assemble1,assemble1) = H11(assemble1,assemble1) + Eps1elem(k,:)'*Eps1elem(k,:);
    end

    % build up A matrix with all constraints
    A = [Aconst; B11];
    % speciality conditions
    B = zeros(size(A,1));

    %% Solving the optimization problem - same as for polynomials
    Q = 1; % initial guesses for stiffness values - required since H relies on Qij

    n_iter = 20; % maximum iterations
    delta_lim = 0.001; % tolerance on optimization

    delta = 10; % starting with a tolerance larger than delta_lim
    i = 1;
    Qold = Q; % stiffnesses from previous iteration required to compute error

    while i<n_iter && delta> delta_lim
        % Hessian matrix
        H = (L*w/n_points)^2*(Q^2)*H11;

        % NOTE: to avoid numerical "Warning: Matrix is close to singular or
        % badly scaled" matrix Opt can be scaled with the parameter corr.
        % It does not change the results of the optimization.
        % To avoid using, put corr = 1
        % corr = max(max(A))/max(max(H)); % normalization coefficient
        corr = 1;
        OptM = [H*corr, A'*corr; A,B]; % matrix for virtual fields optimization

        %vector containing the polynomial coefficients for Q11 and the
        %Lagrange multipliers
        Ya = OptM\Za';

        %remove the lagrangian multipliers from the Y vectors because they are
        %of no interest
        Ya(n_nodes + 1: size(Ya)) = [];

        for k = 1:n_points
            % virtual displacement field
            u1vv = zeros(1,n_nodes);
            assemble1 = assemble(k,:);
            u1vv(assemble1) = u1elem(k,:);
            u1va(k) = u1vv*Ya; % for Q11
        end
        % dynamics code
        %calculating Q11 from the first optimized virtual field
        Q = -material.rho*(L*w/n_points)*(sum(u1va'.*A1));

        % compute difference between the current and previous identified values
        delta = sum((Qold-Q).^2./Q.^2);

        Qold = Q; % store current parameters as old before the next iteration
        i = i+1; % increment the step
    end
    
    %----------------------------------------------------------------------
    % Store diagnostics for the identification
    nIterLog(j) = i;
    mMatrixCond(j) = cond(OptM);
    H = (L*w/n_points)^2*(Q(1)^2)*H11; % used to compute noise sensitivity parameter
    eta11t(j) = sqrt(Ya'*H*Ya); %sensitivity to noise parameter (eta11) for Q11 - coefficient of variation when divided by Q11 (in output statement)
  
    %----------------------------------------------------------------------
    % Store the stiffness found for this frame
    Q11t(j) = Q;
    
    % Reset for the next frame
    B11 = zeros(1,n_nodes);
    u1va = zeros(1,n_points);
    
    %----------------------------------------------------------------------
    % Store the virtual fields for post processing and visualisation
    u1va = zeros(1,n_points);
    Eps1va = zeros(1,n_points);
    for k = 1:n_points
        % virtual displacement field
        u1vv = zeros(1,n_nodes);
        assemble1 = assemble(k,:);
        u1vv(assemble1) = u1elem(k,:);
        u1va(k) = u1vv*Ya; % for Q11

        % Virtual strain fields, 1 components
        Eps1vv=zeros(1,n_nodes);
        Eps1vv(assemble1)=Eps1elem(k,:);
        Eps1va(k)=Eps1vv*Ya; % for Q11
    end
    VFs.u1Star(:,:,j)=reshape(u1va',rows,cols);
    VFs.eps1Star(:,:,j)=reshape(Eps1va',rows,cols);
    u1va = zeros(1,n_points);
    Eps1va = zeros(1,n_points);
    
end

% Push the diagnostics into a data struct to return
VFOptDiag.nIters = nIterLog;
VFOptDiag.mMatrixCond = mMatrixCond; 
VFOptDiag.eta11VsT = eta11t;

% Push the identified stiffnesses into a struct to return
% Q11 approx E11 for this case
identStiffVFOpt.QxxVsT = Q11t;
identStiffVFOpt.ExxVsT = Q11t;

% Turn warnings back on for ill conditioned matrices
warning('on', 'MATLAB:singularMatrix'); 
warning('on', 'MATLAB:nearlySingularMatrix');

end

