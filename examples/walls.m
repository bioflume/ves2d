classdef walls < handle
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
X;          % positions of vesicles
nv;
u;          % velocity field of vesicles
xt;         % tangent unit vector
normal;
sa;         % Jacobian
isa;        % inverse of Jacobian
cur;        % curvature
center;     % center of the point required for stokeslets
            % and rotlets
IK;         % index of Fourier modes for fft and ifft
            % that are needed repetatively
viscCont    % viscosity contrast

end %properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = walls(X,center,u)

% This is the constructor
o.N = size(X,1)/2;              % points per vesicle
o.X = X;                        % position of vesicle
o.nv = 1;
oc = curve;
[o.sa,o.xt,o.cur] = oc.diffProp(o.X);
o.isa = 1./o.sa;
% Jacobian, tangent, and curvature
o.u = u;                % Velocity of vesicle
o.normal = [o.xt(o.N+1:2*o.N,:);-o.xt(1:o.N,:)];
o.center = center;
o.viscCont = zeros(1,1);
% minimum arclength needed for near-singular integration
o.IK = fft1.modes(o.N,1);
% ordering of the fourier modes.  It is faster to compute once here and
% pass it around to the fft differentitation routine

end % wall
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % methods

end % classdef