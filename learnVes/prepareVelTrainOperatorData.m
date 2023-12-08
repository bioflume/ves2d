function prepareVelTrainOperatorData(iset)
%load /workspace/gokberk/relax1step/only1timeStepPCAData.mat
% clear sum of the data to have some space
%clear XnewCoeffs; clear XnewRec; clear XnewStandStore;
%clear XoldCoeffs; clear Xrec; clear colMeans;
%clear errInNew; clear evects;

% load /workspace/gokberk/relax1step/n256Dt1E4Kb1RelaxAllDataSet.mat
load /work2/03353/gokberk/frontera/n256Dt0.0016RelaxAllDataSet.mat
% load ./output/relaxData/completeData/n256Dt0.0016RelaxAllDataSet.mat
clear XnewStandStore;

addpath ../src/
% num. points
N = 128;
op = poten(N);

% iset = 4
nSamples = ones(4,1)*floor(nInstances/4);
nSamples(end) = nInstances-3*nSamples(1);

velMats = zeros(2*N,2*N,nSamples(iset));

idx = 1;
for ives = sum(nSamples(1:iset-1))+1:sum(nSamples(1:iset))
  disp(['Vesicle #' num2str(idx) ' out of ' num2str(nSamples(iset)) ' being processed...'])
  tstart = tic;

  % Build vesicle
  vesicle = capsules(XstandStore(:,ives),[],[],1,1,1);
  vesicle.setUpRate();
  
  % derivatives
  [Ben,Ten,Div] = vesicle.computeDerivs;
  
  % SLP
  G = op.stokesSLmatrix(vesicle);
  
  % Build M and Multiply with Basis
  M = G*Ten*((Div*G*Ten)\eye(vesicle.N))*Div;
  M11 = M(1:end/2,1:end/2); M12 = M(1:end/2,end/2+1:end);
  M21 = M(end/2+1:end,1:end/2); M22 = M(end/2+1:end,end/2+1:end);
  VelMats(:,:,idx) = [M11 M12; M21 M22];
  
  idx = idx + 1;
  tend = toc(tstart);
  disp(['took ' num2str(tend) ' seconds'])
end


fileName = ['/work2/03353/gokberk/frontera/velocityRuns/VelTrainOperatorData_' num2str(iset) '.mat']; 
% fileName = ['./output/relaxData/completeData/veltrainOperatorData_' num2str(iset) '.mat']; 
nsampInSet = nSamples(iset);
save(fileName,'nInstances','nsampInSet','velMats','N','-v7.3')

end
