function prepareNearFieldStokesletData(iset,npar)
load ./workingOnDataSet/output/advectionNetInputX.mat
%clear XnewStandStore;

addpath ../src/
oc = curve;

% Upsampling maybe necessary
nlayers = 3; 
maxLayerDist = @(h) sqrt(h);

% num. points
N = 128;
op = poten(N);
nmodes = 128;

nves = size(XstandStore, 2);

% build basis 
theta = (0:N-1)'/N*2*pi;
ks = (0:N-1)';
basis = 1/N*exp(1i*theta*ks');
Br = real(basis); Bi = imag(basis);

nSamples = ones(npar,1)*floor(nves/npar);
nSamples(end) = nSamples(end) + nves-sum(nSamples);
nSamplesInSet = sum(nSamples(1:iset)) - (sum(nSamples(1:iset-1))+1) + 1;


idx = 1;
for ives = sum(nSamples(1:iset-1))+1:sum(nSamples(1:iset))
  disp(['Vesicle #' num2str(idx) ' out of ' num2str(nSamples(iset)) ' being processed...'])
  tstart = tic;

  % Build vesicle
  vesicle = capsules(XstandStore(:,ives),[],[],1,1,0);
  %vesicle.setUpRate();
  
  % Generate the grid for velocity
  [~,tang] = oc.diffProp(XstandStore(:,ives));
  % get x and y components of normal vector at each point
  nx = tang(N+1:2*N);
  ny = -tang(1:N);

  % Points where velocity is calculated involve the points on vesicle
  tracersX = zeros(2*N, nlayers-1);
  tracersX(:,1) = XstandStore(:,ives);

  % Generate tracers
  h = vesicle.length/vesicle.N;  % arc-length spacing
  dlayer = (0:nlayers-1)'/(nlayers-1) * maxLayerDist(h);
  for il = 2 : nlayers
    tracersX(:,il) = [vesicle.X(1:end/2)+nx*dlayer(il);vesicle.X(end/2+1:end)+ny*dlayer(il)];
  end

  tracers.N = N;
  tracers.nv = nlayers-1;
  tracers.X = tracersX(:,2:nlayers);

  G = op.stokesSLmatrix(vesicle);
  kernel = @op.exactStokesSL;
  kernelDirect = @op.exactStokesSL;
  SLP = @(X) op.exactStokesSLdiag(vesicle,G,X);
  [~,NearV2T] = vesicle.getZone(tracers,2);
  
  VelOnGridModesReal = zeros(2*N,nlayers-1,nmodes);
  VelOnGridModesImag = zeros(2*N,nlayers-1,nmodes);
  selfVelModesReal = zeros(2*N,nmodes);
  selfVelModesImag = zeros(2*N,nmodes);

  for imode = 1 : nmodes
    forRealVels = [Br(:,imode); Bi(:,imode)];
    forImagVels = [-Bi(:,imode); Br(:,imode)];

    VelOnGridModesReal(:,:,imode) = op.nearSingInt(vesicle,forRealVels,SLP,[],NearV2T,kernel,kernelDirect,tracers,false,false);
    selfVelModesReal(:,imode) = G*forRealVels;
    

    VelOnGridModesImag(:,:,imode) = op.nearSingInt(vesicle,forImagVels,SLP,[],NearV2T,kernel,kernelDirect,tracers,false,false);
    selfVelModesImag(:,imode) = G*forImagVels;
  end
  fileName = ['E:\nearFieldData\vesicleByvesicle\vesicleID_' num2str(ives) '.mat']; 
  save(fileName,'VelOnGridModesImag','VelOnGridModesReal','selfVelModesImag','selfVelModesReal','-v7.3')
  
  idx = idx + 1;
  tend = toc(tstart);
  disp(['took ' num2str(tend) ' seconds'])
end


% fileName = ['/work2/03353/gokberk/frontera/velocityRuns/elocityTrain128modesFFTData_' num2str(iset) '.mat']; 
%fileName = ['E:\advectData\veltrain128modesFFTData_' num2str(iset) '.mat']; 
%nsampInSet = nSamples(iset);
%save(fileName,'nInstances','nsampInSet','zRealStore','zImagStore',...
%  'activeModes','N','nmodes','-v7.3')

end
