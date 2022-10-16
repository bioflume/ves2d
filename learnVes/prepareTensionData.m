clear; clc;
%load /workspace/gokberk/relax1step/N96DataSets/n96Relax100KAllDataSet.mat
%load /workspace/gokberk/relax1step/only1timeStepPCAData.mat
load /workspace/gokberk/relax1step/n256Dt1E4Kb1RelaxAllDataSet.mat
clear XnewStandStore;

% clear sum of the data to have some space

%clear XnewCoeffs; clear XnewRec; clear XnewStandStore;
%clear XoldCoeffs; clear Xrec; clear colMeans;
%clear errInNew; clear evects;

addpath ../src/
oc = curve;
% num. points
N = 256;
op = poten(N);
nmodes = 24;

% compute inv(Div*G*Ten)*Div*basis(for 12 vectors), we cannot store all
% compute inv(Div*G*Ten)*Div*G*(-Ben)*x, store 12 modes
% We have standardized shapes in:
% XstandStore = zeros(2*N,nInstances);
% XnewStandStore = zeros(2*N,nInstances);

% build basis 
theta = (0:N-1)'/N*2*pi;
ks = (0:N-1)';
basis = 1/N*exp(1i*theta*ks');
% How to reduce dimension
% depending on activeModes, we choose basis vectors
activeModes = [(1:nmodes/2)';(N-nmodes/2+1:N)'];
% we cannot look at statistics of velocity's energy spectrum
% because we consider arbitrary velocity

iset = 4
nSamples = ones(4,1)*floor(nInstances/4);
nSamples(end) = nInstances-3*nSamples(1);

zRealStore = zeros(N,nmodes,nSamples(iset));
zImagStore = zeros(N,nmodes,nSamples(iset));

benRealStore = zeros(N,nSamples(iset));
benImagStore = zeros(N,nSamples(iset));

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
  
  % Build M and Multiply with Basis (predict Z11 on 48 points)
  M = ((Div*G*Ten)\eye(vesicle.N))*Div;
  M11 = M(:,1:end/2); M12 = M(:,end/2+1:end);
  B1 = real(basis(:,activeModes)); B2 = imag(basis(:,activeModes));
  Z11 = M11*B1+M12*B2; Z12 = M12*B1-M11*B2;
  zRealStore(:,:,idx) = Z11; 
  zImagStore(:,:,idx) = Z12;
   
  if 1
  % Compute the action of M on G*(-Ben)
  benTerm = M*G*(-Ben)*XstandStore(:,ives); 
  benH = fft(benTerm)/N;
  benRealStore(:,idx) = real(benH); % predict 48 modes
  benImagStore(:,idx) = imag(benH);
  end

  % RECONSTRUCTION WORKS AS FOLLOWS
  % PREDICT zReal;zImag for every mode on reduced dimension, i.e.,
  % Z11r = interpft(Z11,24);... So we predict columns of [Z11r Z12r; Z21r
  % Z22r], then parse this matrix as Z11r, Z12r ...
  % V1 = real(coeffs(activeModes)); V2 = imag(coeffs(activeModes));
  % MVinfRed = [Z11r*V1+Z12r*V2;Z21r*V1+Z22r*V2];
  % MVinfFull =
  % [interpft(MVinfRed(1:end/2),N);interpft(MVinfRed(end/2+1:end),N)];
  
  % OLD VERSION (WRONG) Get z = [z1; z2; ...; zk; ...; zn] = M*wk(basis-vectors)
  %   Mz = M(:,1:end/2)+1i*M(:,end/2+1:end);
  %   mult = Mz*basis(:,activeModes);
  %   zRealStore(:,:,idx) = real(mult);
  %   zImagStore(:,:,idx) = imag(mult);
  idx = idx + 1;
  tend = toc(tstart);
  disp(['took ' num2str(tend) ' seconds'])
end


fileName = ['/workspace/gokberk/relax1step/NEWn256trainTenBasis24modesData_' num2str(iset) '.mat']; 
nsampInSet = nSamples(iset);
save(fileName,'nInstances','nsampInSet','zRealStore','zImagStore',...
  'activeModes','N','nmodes','-v7.3')
if 1
fileName = ['/workspace/gokberk/relax1step/NEWn256trainTenBendingData_' num2str(iset) '.mat']; 
nsampInSet = nSamples(iset);
save(fileName,'nInstances','nsampInSet',...
  'activeModes','N','nmodes','benRealStore','benImagStore','-v7.3')
end



