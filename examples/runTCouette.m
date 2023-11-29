VCs =[1; 10];
vSize = 6;
kappa = 1;
RA = 0.6;
vinf = -40;
Xg = 10; % 3+3

parfor i = 1 : 2
VC = VCs(i);
runName = ['./output/size' num2str(vSize) '_Vf' num2str(vinf) '_VC' num2str(VC) '_RA' num2str(RA) '_run'];
driver_rotational_vortex(runName, VC, RA, kappa, vinf, vSize, Xg, pi/2); 
end

