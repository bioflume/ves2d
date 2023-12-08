clear;
addpath '/Users/gokberk/Documents/GitHub/ves2d/learnVes/output/relaxData/completeData'

%load MtimesFFTBasisN96Data_1.mat
%load RA75velocityTrainFFTData_1.mat
% load NEWvelocityTrain24modesFFTData_1.mat
load veltrain128modesFFTData_1.mat
nsampInSets(1) = nsampInSet;
zRealStore1 = zRealStore;
zImagStore1 = zImagStore;
%shapeModesStore1 = shapeModesStore;

%load MtimesFFTBasisN96Data_2.mat
%load RA75velocityTrainFFTData_2.mat
% load NEWvelocityTrain24modesFFTData_2.mat
load veltrain128modesFFTData_2.mat
nsampInSets(2) = nsampInSet;
zRealStore2 = zRealStore;
zImagStore2 = zImagStore;
%shapeModesStore2 = shapeModesStore;

%load MtimesFFTBasisN96Data_3.mat
%load RA75velocityTrainFFTData_3.mat
% load NEWvelocityTrain24modesFFTData_3.mat
load veltrain128modesFFTData_3.mat
nsampInSets(3) = nsampInSet;
zRealStore3 = zRealStore;
zImagStore3 = zImagStore;
%shapeModesStore3 = shapeModesStore;

%load MtimesFFTBasisN96Data_4.mat
%load RA75velocityTrainFFTData_4.mat
% load NEWvelocityTrain24modesFFTData_4.mat
load veltrain128modesFFTData_4.mat
nsampInSets(4) = nsampInSet;
zRealStore4 = zRealStore;
zImagStore4 = zImagStore;
%shapeModesStore4 = shapeModesStore;

%shapeModesStore = zeros(N,2,nInstances);
zRealStore = zeros(2*N,nmodes,nInstances);
zImagStore = zeros(2*N,nmodes,nInstances);

part1s = 1; part1e = nsampInSets(1);
part2s = nsampInSets(1)+1; part2e = sum(nsampInSets(1:2));
part3s = part2e+1; part3e = sum(nsampInSets(1:3));
part4s = part3e+1; part4e = nInstances;

%shapeModesStore(:,:,part1s:part1e) = shapeModesStore1;
%shapeModesStore(:,:,part2s:part2e) = shapeModesStore2;
%shapeModesStore(:,:,part3s:part3e) = shapeModesStore3;
%shapeModesStore(:,:,part4s:part4e) = shapeModesStore4;

zRealStore(:,:,part1s:part1e) = zRealStore1;
zRealStore(:,:,part2s:part2e) = zRealStore2;
zRealStore(:,:,part3s:part3e) = zRealStore3;
zRealStore(:,:,part4s:part4e) = zRealStore4;

zImagStore(:,:,part1s:part1e) = zImagStore1;
zImagStore(:,:,part2s:part2e) = zImagStore2;
zImagStore(:,:,part3s:part3e) = zImagStore3;
zImagStore(:,:,part4s:part4e) = zImagStore4;

save('VelOpFFTBasisN128All128modesData.mat','nInstances',...
  'zRealStore','zImagStore','activeModes','N','nmodes','-v7.3')
