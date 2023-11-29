function runManyTGs
VCs = [4;5];
vSize = 6;
kappa = 1;
RA = 0.6;
VFs = [400; 800; 1600; 2000];

[vv, vf] = meshgrid(VCs,VFs);
parpool(8);
parfor irun = 1 : numel(vv(:))
  VC = vv(irun);
  VF = vf(irun);
  runName = ['/work2/03353/gokberk/frontera/tgRuns/Vsize' num2str(vSize) '_Vf' num2str(VF) '_VC' num2str(VC) '_RA' num2str(RA)  '_run'];
  driver_taylor_green(runName, VC, RA, kappa, VF, vSize, 0.25, pi/2); 
end
end

