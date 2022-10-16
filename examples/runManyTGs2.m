function runManyTGs2(RA,vSize,Xg,IA,iset)
kappa = 1;
%vcs = [5;10;20;50;75;100];
%vinfs = [8;32;64;128;256];

%[VC, VF] = meshgrid(vcs,vinfs);
VC = [4; 5; 6; 7; 2; 3; 4];
VF = [192; 384; 384; 384; 0.25; 0.25; 0.25];

parpool(7);
parfor irun = 1 : numel(VF(:))
  runName = ['/mnt/ceph/users/gkabacaoglu/SVTGRuns/set14/Vsize' num2str(vSize) '_Vf' num2str(VF(irun)) '_VC' num2str(VC(irun)) '_RA' num2str(RA) '_Xg' num2str(Xg) '_IA' num2str(IA) '_run'];
  driver_taylor_green(runName, VC(irun), RA, kappa, VF(irun), vSize, Xg, IA); 
end
end

