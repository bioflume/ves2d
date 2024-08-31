start_mode = [1; 33; 65; 97];
end_mode = [32; 64; 96; 128];
clear; clc;
load(['TenBasis128modesData_1'])
nfiles = 4;
zRealStore1_32 = zeros(128,32,nInstances);
zImagStore1_32 = zeros(128,32,nInstances);
nLast = 0;
for ff = 1 : nfiles
  load(['TenBasis128modesData_' num2str(ff)])
  
  if ff == 1; save('inputToAdvectTensionNets.mat','XstandStore'); end;
  clear XstandStore
  
  zRealStore1_32(:,:,nLast+1:nLast+nsampInSet) = zRealStore(:,1:32,:); % 1st mode is zero
  zImagStore1_32(:,:,nLast+1:nLast+nsampInSet) = zImagStore(:,1:32,:); % 1st mode is zero
  clear zRealStore zImagStore
  nLast = nLast + nsampInSet;
end
save('TenAdvectFFTBasisN128modes1to32Data.mat', 'zRealStore1_32','zImagStore1_32','nInstances','-v7.3')
clear zRealStore1_32 zImagStore1_32
disp('first set is done')
pause
%%
zRealStore33_64 = zeros(128,32,nInstances);
zImagStore33_64 = zeros(128,32,nInstances);
nLast = 0;
for ff = 1 : nfiles
  load(['TenBasis128modesData_' num2str(ff)])
  
  
  zRealStore33_64(:,:,nLast+1:nLast+nsampInSet) = zRealStore(:,33:64,:); % 1st mode is zero
  zImagStore33_64(:,:,nLast+1:nLast+nsampInSet) = zImagStore(:,33:64,:); % 1st mode is zero
  clear zRealStore zImagStore
  nLast = nLast + nsampInSet;
end
save('TenAdvectFFTBasisN128modes33to64Data.mat', 'zRealStore33_64','zImagStore33_64','nInstances','-v7.3')
clear zRealStore33_64 zImagStore33_64


zRealStore65_96 = zeros(128,32,nInstances);
zImagStore65_96 = zeros(128,32,nInstances);
nLast = 0;
for ff = 1 : nfiles
  load(['TenBasis128modesData_' num2str(ff)])
  
  
  zRealStore65_96(:,:,nLast+1:nLast+nsampInSet) = zRealStore(:,65:96,:); % 1st mode is zero
  zImagStore65_96(:,:,nLast+1:nLast+nsampInSet) = zImagStore(:,65:96,:); % 1st mode is zero
  clear zRealStore zImagStore
  nLast = nLast + nsampInSet;
end
save('TenAdvectFFTBasisN128modes65to96Data.mat', 'zRealStore65_96','zImagStore65_96','nInstances','-v7.3')
clear zRealStore65_96 zImagStore65_96


zRealStore97_128 = zeros(128,32,nInstances);
zImagStore97_128 = zeros(128,32,nInstances);
nLast = 0;
for ff = 1 : nfiles
  load(['TenBasis128modesData_' num2str(ff)])
  
  
  zRealStore97_128(:,:,nLast+1:nLast+nsampInSet) = zRealStore(:,97:128,:); % 1st mode is zero
  zImagStore97_128(:,:,nLast+1:nLast+nsampInSet) = zImagStore(:,97:128,:); % 1st mode is zero
  clear zRealStore zImagStore
  nLast = nLast + nsampInSet;
end
save('TenAdvectFFTBasisN128modes97to128Data.mat', 'zRealStore97_128','zImagStore97_128','nInstances','-v7.3')
clear zRealStore97_128 zImagStore97_128
