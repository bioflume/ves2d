clear; clc;
nfiles = 4;
nInstances = 156225;


for ff = 1 : nfiles
  %load(['./advectData/veltrain128modesFFTData_' num2str(ff)])
  load(['./advectTensionData/TenBasis128modesData_' num2str(ff) '.mat'])
  for imode = 1 : 128
    zRealStorePart = reshape(zRealStore(:,imode,:),128,nsampInSet);
    zImagStorePart = reshape(zImagStore(:,imode,:),128,nsampInSet);
    save(['./combinedTensionAdvectAllModes/partsModes/TensionAdvectNetFFTBasisN128_mode' num2str(imode) '_part' num2str(ff) '_Data.mat'], 'zRealStorePart','zImagStorePart','nsampInSet','-v7.3')
  end
  disp(['File ' num2str(ff) ' is done'])
end
