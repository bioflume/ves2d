clear; 
kappa = 1e-1;
RA = 0.6;
VF = 8;
vSize = 9;
[VCs, Xgs] = meshgrid([1;5;50],[1/16;1/8;1/6]);
parpool(9);
parfor irun = 1 : numel(VCs(:))
  runName = ['/mnt/ceph/users/gkabacaoglu/SVTGRuns/set3/Xg' num2str(Xgs(irun)) '_VC' num2str(VCs(irun)) '_run'];
  driver_taylor_green(runName, VCs(irun), RA, kappa, VF, vSize, Xgs(irun)); 
end


