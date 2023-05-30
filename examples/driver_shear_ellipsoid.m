clear;

N = [16; 32; 64; 128; 256];

ave_iters = [];
errors = [];
for ir = 1 : numel(N)

% [error, ave_iter, orientation, theta_theo] = ellipsoid_shearTestSymmAlpert(1,N(ir)); 
[error, ave_iter, orientation, theta_theo] = ellipsoid_shearTestSDDj(1,N(ir)); 
errors(:,ir) = error;
ave_iters(ir,1) = ave_iter;
final_orientation(ir) = orientation;
% orientation
% theta_theo
% pause
end

final_orientation = pi - final_orientation - pi/2; % theoretical value for short semiaxis
error_in_final_orient = abs(theta_theo + final_orientation)/abs(theta_theo);