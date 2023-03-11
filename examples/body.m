classdef body < handle
% This class implements standard calculations that need to
% be done to a vesicle, solid wall, or a collection of arbitrary
% target points (such as tracers or pressure/stress targets)
% Given a vesicle, the main tasks that can be performed are
% computing the required derivatives (bending, tension, surface
% divergence), the traction jump, the pressure and stress, 
% and constructing structures required for near-singluar
% integration

properties
N;          % number of points per vesicle
nv;         % number of vesicles
X;          % positions of vesicles
traction;        % traction of vesicles
u;          % velocity field of vesicles
radius;     % radius
xt;         % tangent unit vector
normal;
sa;         % Jacobian
isa;        % inverse of Jacobian
length;     % minimum length over all vesicles
cur;        % curvature
center;     % center of the point required for stokeslets
            % and rotlets
IK;         % index of Fourier modes for fft and ifft
            % that are needed repetatively
orientation % orientation of body
viscCont    % viscosity contrast

end %properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = body(X,radius,center,orientation,u)

% This is the constructor
o.N = size(X,1)/2;              % points per vesicle
o.nv = size(X,2);               % number of vesicles
o.X = X;                        % position of vesicle
oc = curve;
[o.sa,o.xt,o.cur] = oc.diffProp(o.X);
o.isa = 1./o.sa;
% Jacobian, tangent, and curvature
o.traction = zeros(size(X));          % Traction of vesicle
o.u = u;                % Velocity of vesicle
o.normal = [o.xt(o.N+1:2*o.N,:);-o.xt(1:o.N,:)];
o.center = center;
o.radius = radius;
o.orientation = orientation;
o.viscCont = zeros(o.nv,1);
[~,~,len] = oc.geomProp(X);
o.length = min(len);
% minimum arclength needed for near-singular integration
o.IK = fft1.modes(o.N,o.nv);
% ordering of the fourier modes.  It is faster to compute once here and
% pass it around to the fft differentitation routine

end % body
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function K = calc_K_matrix(o)
% K matrix is a matrix of size (2*N,3,nv)
% for velocity [u, v, w]; w > 0 for CCW direction

K = zeros(2*o.N, 3,o.nv);

for iv = 1 : o.nv
  cx = mean(o.X(1:o.N,iv));
  cy = mean(o.X(o.N+1:end,iv));  
  for ij = 1 : o.N
    x = o.X(ij,iv)-cx; y = o.X(ij+o.N,iv)-cy;
    K(ij,:,iv) = [1 0 y];
    K(ij+o.N,:,iv) = [0 1 -x];
  end% ij
end% iv

end % calc_K_matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function K = calc_K_matrix_weighted(o)
% K matrix is a matrix of size (2*N,3,nv)
% for velocity [u, v, w]; w > 0 for CCW direction

K = zeros(2*o.N, 3,o.nv);
sa = o.sa; % jacobian

for iv = 1 : o.nv
  cx = mean(o.X(1:o.N,iv));
  cy = mean(o.X(o.N+1:end,iv));  
  for ij = 1 : o.N
    x = o.X(ij,iv)-cx; y = o.X(ij+o.N,iv)-cy;
    K(ij,:,iv) = [1 0 y];
    K(ij+o.N,:,iv) = [0 1 -x];
  end% ij
end% iv

sa = [sa;sa];
K = K .* sa(:,ones(3,1)) * 2 * pi / o.N;

end % calc_K_matrix_weighted

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function K = calc_KT_matrix(o,sa)
% Adjoint of K matrix

K = zeros(3,2*o.N,o.nv);

for iv = 1 : o.nv
  cx = mean(o.X(1:o.N,iv));
  cy = mean(o.X(o.N+1:end,iv));  
  for ij = 1 : o.N
    x = o.X(ij,iv)-cx; y = o.X(ij+o.N,iv)-cy;
    K(:,ij,iv) = [1;0;y];
    K(:,ij+o.N,iv) = [0;1;-x];
  end% ij
end% iv
sa = [sa;sa]';
K = K .* sa(ones(3,1),:) * 2*pi/o.N;

end % calc_KT_matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function K = calc_KT_matrix_weightless(o,sa)
% Adjoint of K matrix

K = zeros(3,2*o.N,o.nv);

for iv = 1 : o.nv
  cx = mean(o.X(1:o.N,iv));
  cy = mean(o.X(o.N+1:end,iv));  
  for ij = 1 : o.N
    x = o.X(ij,iv)-cx; y = o.X(ij+o.N,iv)-cy;
    K(:,ij,iv) = [1;0;y];
    K(:,ij+o.N,iv) = [0;1;-x];
  end% ij
end% iv

end % calc_KT_matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % methods

end % classdef