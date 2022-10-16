function runManyTGs(vSize, VC)
kappa = 1;
ras = [0.3; 0.6; 0.9];
vinfs = [1e-1; 5e-1; 1; 2];


[RA, VF] = meshgrid(ras,vinfs);
parpool(12);
parfor irun = 1 : numel(RA(:))
  runName = ['/mnt/ceph/users/gkabacaoglu/SVTGRuns/set6/Vsize' num2str(vSize) '_Vf' num2str(VF(irun)) '_VC' num2str(VC) '_RA' num2str(RA(irun)) '_run'];
  driver_taylor_green(runName, VC, RA(irun), kappa, VF(irun), vSize, 0.05); 
end
end

