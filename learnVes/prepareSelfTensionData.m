function prepareSelfTensionData(iset,npar)
%load /workspace/gokberk/relax1step/N96DataSets/n96Relax100KAllDataSet.mat
%load /workspace/gokberk/relax1step/only1timeStepPCAData.mat
%load /workspace/gokberk/relax1step/n256Dt1E4Kb1RelaxAllDataSet.mat
%clear XnewStandStore;

load ./necessaryMatFiles/X106KinitShapes.mat
% Xstore Loaded
% clear sum of the data to have some space

%clear XnewCoeffs; clear XnewRec; clear XnewStandStore;
%clear XoldCoeffs; clear Xrec; clear colMeans;
%clear errInNew; clear evects;

addpath ../src/
oc = curve;
% num. points
N = 128;
op = poten(N);

% store aligned shapes
XstandStore = [];
nInstances = size(Xstore,2);
% 


nSamples = ones(npar,1)*floor(nInstances/npar);
nSamples(end) = nInstances-(npar-1)*nSamples(1);
dnn = dnnTools;

selfTenStore = [];

idx = 1;
for ives = sum(nSamples(1:iset-1))+1:sum(nSamples(1:iset))
  disp(['Vesicle #' num2str(idx) ' out of ' num2str(nSamples(iset)) ' being processed...'])
  tstart = tic;
  % change the resolution
  Xinit = [interpft(Xstore(1:end/2,ives),N); interpft(Xstore(end/2+1:end,ives),N)]; 
  
  % standardize
  [Xinit,scaling,rotate,trans,sortIdx] = dnn.standardizationStep(Xinit,oc);

  XstandStore(:,idx) = Xinit;

  % Build vesicle
  vesicle = capsules(Xinit,[],[],1,1,1);
  vesicle.setUpRate();
  
  % derivatives
  [Ben,Ten,Div] = vesicle.computeDerivs;
  
  % SLP
  G = op.stokesSLmatrix(vesicle);
  
  % Build M and Multiply with Basis (predict Z11 on 48 points)
  M = ((Div*G*Ten)\eye(vesicle.N))*Div;
  
  % Compute the action of M on G*(-Ben)
  selfTenStore(:,idx) = M*G*(-Ben)*Xinit; 
  
  idx = idx + 1;
  tend = toc(tstart);
  disp(['took ' num2str(tend) ' seconds'])
end


fileName = ['./output/selfTensionData/selfTensionData_' num2str(iset) '.mat']; 
nsampInSet = nSamples(iset);
save(fileName,'nInstances','nsampInSet','XstandStore','selfTenStore','-v7.3')



