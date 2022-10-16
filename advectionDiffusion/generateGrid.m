function [x,y,r1D,theta,s,R1D,l] = generateGrid(N_theta,N_radii,radius)

% Fourier grid
theta = [0:N_theta-1]'/(N_theta)*2*pi; 

% Chebyshev grid
% s   = [0:N_radii]'*pi/N_radii; %[0,pi]
% R1D = -cos(s); % R1D in [-1,1]
l   = [0:N_radii]';

%% New Chebyshev grid
s = (2*[1 : (N_radii+1)]' - 1) * pi / (2*(N_radii+1)); % (0,pi)
R1D = -cos(s); %(-1,1)
%%
% Linear map to exact geometry
% r1D = (R1D+1)*(radius(2)-radius(1))/2 + radius(1); % for old Cheb grid
r1D = (R1D-R1D(1))*(radius(2)-radius(1))/(R1D(end)-R1D(1)) + radius(1); % r1D in [r1,r2] 

% Mesh
% [THETA,R]  = meshgrid([theta;2*pi],R1D); % in meshgrid form
[THETA,r]      = meshgrid([theta;2*pi],r1D); % [r1,r2] in meshgrid form

x = r.*cos(THETA); y = r.*sin(THETA);          % in cartesian coordinates
end