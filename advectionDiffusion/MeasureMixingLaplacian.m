function [H1Norm] = MeasureMixingLaplacian(C,N_radii,N_theta,radius,BC_ND)

% ************************************
BC_ND = 2; % Apply Dirichlet BC for H-1 norm
% ************************************

% Interpolate for spectral coefficients
[~, ~, ~, ~, a0, ac, as, aM] = interpFFT(N_theta, N_radii, C);

%% ------------------------------------------------------------------------
% SPACE DISCRETIZATION
% -------------------------------------------------------------------------
M = floor((N_theta+1)/2); % number of Fourier coefficients for each (cos and sin) (~k)

[~,~,r1D,theta,s,R1D,l] = generateGrid(N_theta,N_radii,radius);

%% Setting the stage for the calculations
% Trig. functions
bcos  = @(l,s) cos(s*l');
bsin  = @(l,s) sin(s*l');

% Derivatives of "s" w.r.t. Chebyshev points (R)
dsdR = @(R) 1./sqrt(1-R.^2);
d2sdR2 = @(R) R./sqrt((1-R.^2).^3);

% Derivative of R w.r.t. exact geometry r (due to linear mapping)
dRdr = (R1D(end)-R1D(1))/(radius(2)-radius(1));

%% Forming the system of equations 
% Multipliers of Chebyshev Coefficients (b^c_kl, b_ol) (LHS) = (I-Laplacian)
BC_Multip = @(l,s,dsdR,d2sdR2,r,sgn) sgn*(repmat(dsdR.^2,1,length(l)).*bcos(l,s).*repmat(l'.^2,length(r),1)*dRdr^2 ...
                            + repmat(d2sdR2,1,length(l)).*bsin(l,s).*repmat(l',length(r),1)*dRdr^2 ...
                            + repmat(dsdR,1,length(l)).*repmat(1./r,1,length(l)).*bsin(l,s).*repmat(l',length(r),1)*dRdr) ...
                            + bcos(l,s);

% Compute LHS matrices (excluding BC nodes, i.e. r1 and r2) for Cheb. pts.
% Therefore, constant sections of the matrices are computed once
% (LHS_S = LHS_C, thus, is not created)
LHS_C = zeros(N_radii+1,N_radii+1);

LHS_C(2:end-1,:) = BC_Multip(l,s(2:end-1),dsdR(R1D(2:end-1)),d2sdR2(R1D(2:end-1)),r1D(2:end-1), 1);

% Boundary Nodes
if BC_ND == 1
    % New Cheb Grid (Neumann BCs)
    LHS_C(1,:) = (sin(s(1)*l)).*l;
    LHS_C(N_radii+1,:) = (sin(s(end)*l)).*l;
else
    % Dirichlet BCs
    LHS_C(1,:) = cos(s(1)*l);
    LHS_C(end,:) = cos(s(end)*l);
end

w1 = 1/sqrt(N_radii+1);
wL = sqrt(2/(N_radii+1));

LHS_C(:,1) = w1 * LHS_C(:,1); 
LHS_C(:,2:end) = wL * LHS_C(:,2:end); 

% Fourier Mode Matrix
FourModeMatrix = repmat(1./(r1D(2:end-1).^2),1,length(l)).*bcos(l,s(2:end-1)); 
FourModeMatrix(:,1) = w1 * FourModeMatrix(:,1); 
FourModeMatrix(:,2:end) = wL * FourModeMatrix(:,2:end);

% Initialize the Chebyshev Coeff Matrices for COS and SIN
B_C = zeros(N_radii+1,M-1);
B_S = zeros(N_radii+1,M-1);

for kk = 1 : M-1
    
    % For ac and as:
    %----------------------------------------------------------------------
    % Form LHS:
    LHS_CK = LHS_C;
    LHS_CK(2:end-1,:) = LHS_CK(2:end-1,:) + kk^2*FourModeMatrix;
    
    % Form RHS:
    RHS_FC = ac(:,kk);
    RHS_FS = as(:,kk);

    % Apply BCs:
    RHS_FC([1;N_radii+1]) = [0;0];
    RHS_FS([1;N_radii+1]) = [0;0];

    % Solve:
    B_C(:,kk) = LHS_CK\RHS_FC;
    B_S(:,kk) = LHS_CK\RHS_FS;

end

% For am:
%--------------------------------------------------------------------------
% Form LHS:
LHS_CM = LHS_C;
LHS_CM(2:end-1,:) = LHS_CM(2:end-1,:) + M^2*FourModeMatrix;

% Form RHS:
RHS_FM = aM;

% Apply BCs:
RHS_FM([1;N_radii+1]) = [0;0];  

% Solve:
B_M = LHS_CM\RHS_FM;

% For a0:
%--------------------------------------------------------------------------
% Form RHS:
RHS_F0 = a0;

% Apply BCs:
RHS_F0([1;N_radii+1]) = [0;0];

% Solve:
B_0 = LHS_C\RHS_F0;

% Reconstruct Concentration
G = reconst_fun(N_theta, N_radii, B_0, B_C, B_S, B_M); 
% CG = C.*G;

rMat = repmat(r1D,1,size(G,2));

H1Norm = sqrt(trapz(theta, trapz(r1D, C.*G.*rMat,1)));

