clear; clc;

shearStrength = 1;
dt = 0.05;
Ns = [16;32;64;128;256];
% [ave_iterA,max_iterA,trajectoriesA] = starSuspension_TestSymmAlpert(shearStrength, dt, Ns(3));
[ave_iter,max_iter,trajectories] = starSuspension_Test1stKind(shearStrength, dt, Ns(3));

% [ave_iter,max_iter,trajectories] = starSuspensionConfined_TestSymmAlpert(shearStrength, dt, Ns(3));

% [ave_iterP,max_iterP,trajectoriesP] = starSuspension_TestPowerMiranda(shearStrength, dt, Ns(3));