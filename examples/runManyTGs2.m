function runManyTGs2(RA,vSize,Xg,IA)
kappa = 1;
vcs = [1;2;3;4;5;8;10;15];
vinfs = [0.1;0.5;2;8];
%vinfs = [16; 32; 48; 64];
[VC, VF] = meshgrid(vcs,vinfs);
parpool(32);
parfor irun = 1 : numel(VF(:))
  runName = ['/mnt/ceph/users/gkabacaoglu/SVTGRuns/set8/Vsize' num2str(vSize) '_Vf' num2str(VF(irun)) '_VC' num2str(VC(irun)) '_RA' num2str(RA) '_Xg' num2str(Xg) '_IA' num2str(IA) '_run'];
  driver_taylor_green(runName, VC(irun), RA, kappa, VF(irun), vSize, Xg, IA); 
end
end

