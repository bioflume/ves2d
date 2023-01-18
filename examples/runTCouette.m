VC = 1;
vSize = 6;
kappa = 1;
RA = 0.7;
vinf = -10;
Xg = 10; % 3+3
runName = ['./output/size' num2str(vSize) '_Vf' num2str(vinf) '_VC' num2str(VC) '_RA' num2str(RA) '_run'];
driver_rotational_vortex(runName, VC, RA, kappa, vinf, vSize, Xg, pi/2); 


