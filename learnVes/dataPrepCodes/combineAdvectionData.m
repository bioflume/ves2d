clear; clc;
nfiles = 4;
nInstances = 156225;

for imode = 1 : 128
  zRealStore = zeros(128,nInstances);
  zImagStore = zeros(128,nInstances);
  nLast = 0;
  for ff = 1 : nfiles
    load(['./combinedTensionAdvectAllModes/partsModes/TensionAdvectNetFFTBasisN128_mode' num2str(imode) '_part' num2str(ff) '_Data.mat'])
    zRealStore(:,nLast+1:nLast+nsampInSet) = zRealStorePart;%reshape(zRealStore(:,imode,:),256,nsampInSet); % 1st mode is zero
    zImagStore(:,nLast+1:nLast+nsampInSet) = zImagStorePart;%reshape(zImagStore(:,imode,:),256,nsampInSet); % 1st mode is zero
    nLast = nLast + nsampInSet;
  end
  save(['./combinedTensionAdvectAllModes/TensionAdvectNetFFTBasisN128_mode' num2str(imode) 'Data.mat'], 'zRealStore','zImagStore','nInstances','-v7.3')
  disp(['Mode ' num2str(imode) ' is done'])
end

