clear; clc;
addpath ../src/

oc = curve;
dnn = dnnTools;

% if shapes are already standardized
istandard = true;

% load shapes
%load /workspace/gokberk/X100KinitShapes.mat
%load /workspace/gokberk/relax1step/n96Relax100KAllDataSet.mat
%load /workspace/gokberk/relax1step/n256FineDtRelax100KDataSet.mat
%load /workspace/gokberk/relax1step/only1timeStepPCAData.mat
load /workspace/gokberk/relax1step/n256Dt2p5E3Kb1RelaxAllDataSet.mat 
kappa = 1;

% XstandStore, XnewStandStore

% num. points
N = 256;

% number of principal components
nComp = 64;

% if not standardized, then do so
if ~istandard
Xstore = XstandStore;
XnewStore = XnewStandStore;
XstandStore = zeros(2*N,nInstances);
XnewStandStore = zeros(2*N,nInstances);

for k = 1:nInstances
  disp([num2str(k) 'th vesicle is running'])
  
  [XstandStore(:,k),scaling,rotate,trans,sortIdx] = ...
    dnn.standardizationStep(Xstore(:,k),oc);

  % Assuming that Xstore(:,k) has equal distribution of points along arc-length
  % so standardize, XnewStore based on Xstore
  XnewStandStore(:,k) = ...
    dnn.standardize(XnewStore(:,k),trans,rotate,scaling,sortIdx,256);
end
end % ~istandard

% PCA
%Xmat = XstandStore';
%colMeans = mean(Xmat);
%[evects,XoldCoeffs] = pca(Xmat);
load ./necessaryMatFiles/pcaBasisNewest


% get coeffs for new shapes
XoldCoeffs = zeros(nInstances,nComp);
XnewCoeffs = zeros(nInstances,nComp);
for k = 1 : nInstances
  XoldCoeffs(k,:) = (XstandStore(:,k)'-colMeans)*evects(:,1:nComp);
  XnewCoeffs(k,:) = (XnewStandStore(:,k)'-colMeans)*evects(:,1:nComp);
end
XnewRec = XnewCoeffs*evects(:,1:nComp)';
XnewRec = bsxfun(@plus, XnewRec, colMeans)';

XoldRec = XoldCoeffs(:,1:nComp)*evects(:,1:nComp)';
XoldRec = bsxfun(@plus, XoldRec, colMeans)';

errInOld = zeros(nInstances,1);
errInNew = zeros(nInstances,1);
for k = 1 : nInstances
  errInOld(k) = sqrt(1/N*sum((XoldRec(1:end/2,k)-XstandStore(1:end/2,k)).^2+...
      (XoldRec(end/2+1:end,k)-XstandStore(end/2+1:end,k)).^2))/sqrt(...
      1/N*sum(XstandStore(:,k).^2));   

  errInNew(k) = sqrt(1/N*sum((XnewRec(1:end/2,k)-XnewStandStore(1:end/2,k)).^2+...
      (XnewRec(end/2+1:end,k)-XnewStandStore(end/2+1:end,k)).^2))/sqrt(...
      1/N*sum(XnewStandStore(:,k).^2));    
end

fileName = '/workspace/gokberk/relax1step/n256Dt2p5E3Kb1nModes64PCAData.mat'; 
save(fileName,'nInstances','colMeans','evects','XnewCoeffs','XoldCoeffs','kappa','dt')

%save('pcaBasisNewest.mat','colMeans','evects')
%fileName = '/workspace/gokberk/relax1step/only1timeStepN96Dt_PCAerror.mat';
%save(fileName,'nInstances','XnewRec','XoldRec','XstandStore','XnewStandStore',...
%'errInNew','errInOld')

