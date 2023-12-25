if 1
% For the first 32 modes load:
load TenAdvectFFTBasisN128modes1to32Data.mat
end
if 0
% For the modes between 33 and 64 load:
load TenAdvectFFTBasisN128modes33to64Data.mat
end
if 0
% For the modes between 65 and 96 load:
load TenAdvectFFTBasisN128modes65to96Data.mat
end
if 0
% For the modes between 97 and 128 load:
load TenAdvectFFTBasisN128modes97to128Data.mat
end

% these folders include 
% zRealStore1_32 (128, 32 modes, nSamples)
% zImagStore1_32 (128, 32 modes, nSamples)
% for every 32 modes
% there must be ~100K samples in these sets


% The input is stored in XstandStore in the following file
load inputToAdvectTensionNets.mat.mat XstandStore


%  So the networks approximate the following maps
% for every imode = 1 to 128
% XstandStore--> zRealStore(:,imode,nInstances); zImagStore(:,imode,nInstances)



