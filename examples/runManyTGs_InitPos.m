clear; 
kappa = 1;
RA = 0.6;
VF = 128;
vSize = 6;
VC = 10;
Xgs = [0.4];
%Xgs = [0.005; 0.01; 0.02;0.05;0.1;0.15;0.2;0.25;0.3;0.35;0.4];
%Xg = [0.1];
%Xgs = [0.3;0.4];
%IAs = [0;pi/8;pi/6;pi/4;pi/3;pi/2;2*pi/3;3*pi/4;5*pi/6;7*pi/8];
IA = pi/2;
parpool(1);
parfor irun = 1 : numel(Xgs(:))
  runName = ['/mnt/ceph/users/gkabacaoglu/SVTGRuns/set11/LONGXg' num2str(Xgs(irun)) '_VC' num2str(VC) '_run'];
  driver_taylor_green(runName, VC, RA, kappa, VF, vSize, Xgs(irun),IA); 
end


