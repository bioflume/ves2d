clear;
if 0
%load n256trainTenBasisData_1.mat
load NEWn256trainTenBasis24modesData_1.mat
nsampInSets(1) = nsampInSet;
zRealStore1 = zRealStore;
zImagStore1 = zImagStore;
%shapeModesStore1 = shapeModesStore;

%load n256trainTenBasisData_2.mat
load NEWn256trainTenBasis24modesData_2.mat
nsampInSets(2) = nsampInSet;
zRealStore2 = zRealStore;
zImagStore2 = zImagStore;
%shapeModesStore2 = shapeModesStore;

%load n256trainTenBasisData_3.mat
load NEWn256trainTenBasis24modesData_3.mat
nsampInSets(3) = nsampInSet;
zRealStore3 = zRealStore;
zImagStore3 = zImagStore;
%shapeModesStore3 = shapeModesStore;

%load n256trainTenBasisData_4.mat
load NEWn256trainTenBasis24modesData_4.mat
nsampInSets(4) = nsampInSet;
zRealStore4 = zRealStore;
zImagStore4 = zImagStore;
%shapeModesStore4 = shapeModesStore;

zRealStore = zeros(N,nmodes,nInstances);
zImagStore = zeros(N,nmodes,nInstances);

part1s = 1; part1e = nsampInSets(1);
part2s = nsampInSets(1)+1; part2e = sum(nsampInSets(1:2));
part3s = part2e+1; part3e = sum(nsampInSets(1:3));
part4s = part3e+1; part4e = nInstances;

zRealStore(:,:,part1s:part1e) = zRealStore1;
zRealStore(:,:,part2s:part2e) = zRealStore2;
zRealStore(:,:,part3s:part3e) = zRealStore3;
zRealStore(:,:,part4s:part4e) = zRealStore4;

zImagStore(:,:,part1s:part1e) = zImagStore1;
zImagStore(:,:,part2s:part2e) = zImagStore2;
zImagStore(:,:,part3s:part3e) = zImagStore3;
zImagStore(:,:,part4s:part4e) = zImagStore4;

save('NEW_tenMatTimesFFTBasisN256All24modesData.mat','nInstances',...
  'zRealStore','zImagStore','activeModes','N','nmodes','-v7.3')
end

if 1
disp('Combining tension on bending data')
clear;

load NEWn256trainTenBendingData_1.mat
nsampInSets(1) = nsampInSet;
zRealStore1 = benRealStore;
zImagStore1 = benImagStore;

load NEWn256trainTenBendingData_2.mat
nsampInSets(2) = nsampInSet;
zRealStore2 = benRealStore;
zImagStore2 = benImagStore;

load NEWn256trainTenBendingData_3.mat
nsampInSets(3) = nsampInSet;
zRealStore3 = benRealStore;
zImagStore3 = benImagStore;

load NEWn256trainTenBendingData_4.mat
nsampInSets(4) = nsampInSet;
zRealStore4 = benRealStore;
zImagStore4 = benImagStore;

benRealStore = zeros(N,nInstances);
benImagStore = zeros(N,nInstances);

part1s = 1; part1e = nsampInSets(1);
part2s = nsampInSets(1)+1; part2e = sum(nsampInSets(1:2));
part3s = part2e+1; part3e = sum(nsampInSets(1:3));
part4s = part3e+1; part4e = nInstances;

benRealStore(:,part1s:part1e) = zRealStore1;
benRealStore(:,part2s:part2e) = zRealStore2;
benRealStore(:,part3s:part3e) = zRealStore3;
benRealStore(:,part4s:part4e) = zRealStore4;

benImagStore(:,part1s:part1e) = zImagStore1;
benImagStore(:,part2s:part2e) = zImagStore2;
benImagStore(:,part3s:part3e) = zImagStore3;
benImagStore(:,part4s:part4e) = zImagStore4;

save('NEW_tenMatOnBendN256AllData.mat','nInstances',...
  'benRealStore','benImagStore','activeModes','N','nmodes','-v7.3')
end
