function prepareNearFieldVelocityData(iset,npar)

% Load Xstore containing nves number of vesicles 
% Xstore holds 2*N entries in its columns as [x;y] coordinates of nves
% number of vesicles
load ./necessaryMatFiles/X106KinitShapes.mat 
%/work2/03353/gokberk/frontera/X100KinitShapes.mat
addpath ../src/

oc = curve;

% Information about loaded data
nves = size(Xstore,2); % number of vesicles in the data set
Nstore = size(Xstore,1)/2; % number of discretization points on vesicles

% Upsampling maybe necessary
Nup = 128;
nlayers = 5; 
maxLayerDist = @(h) sqrt(h);

% store aligned shapes -- we standardize vesicle shapes
XstandStore = [];


% poten operator
op = poten(Nup);
dnn = dnnTools;

% Procedure: 
% 1. Given the shape, calculate tension and forces
% 2. then calculate velocity at target points


% Divide the data into npar data sets
nSamples = ones(npar,1)*floor(nves/npar);
nSamples(end) = nSamples(end) + nves-sum(nSamples);
nSamplesInSet = sum(nSamples(1:iset)) - (sum(nSamples(1:iset-1))+1) + 1;
nearVelocity = zeros(2*Nup,nlayers,nSamplesInSet);
tracersXstore = zeros(2*Nup,nlayers,nSamplesInSet);
XstandStore = zeros(2*Nup,nSamplesInSet);

idx = 1;
for k = sum(nSamples(1:iset-1))+1 : sum(nSamples(1:iset))
  disp(['Vesicle #' num2str(idx) ' out of ' num2str(nSamples(iset)) ' being processed...'])
  tstart = tic;

  % change the resolution
  Xinit = [interpft(Xstore(1:end/2,k),Nup); interpft(Xstore(end/2+1:end,k),Nup)]; 
  % Standardize if not done yet
  [Xinit,scaling,rotate,trans,sortIdx] = dnn.standardizationStep(Xinit,oc);
  
  [~,area,len] = oc.geomProp(Xinit);
  
  % store the shape
  XstandStore(:,idx) = Xinit;
  
  % Generate the grid for velocity
  [~,tang] = oc.diffProp(Xinit);
  % get x and y components of normal vector at each point
  nx = tang(Nup+1:2*Nup);
  ny = -tang(1:Nup);

  % Points where velocity is calculated involve the points on vesicle
  tracersX = zeros(2*Nup, nlayers-1);
  tracersX(:,1) = Xinit;

  % initialize vesicle
  vesicle = capsules(Xinit, [], [], 1, 1, 1);
  vesicle.setUpRate();

  % Generate tracers
  h = vesicle.length/vesicle.N;  % arc-length spacing
  dlayer = (0:nlayers-1)'/(nlayers-1) * maxLayerDist(h);
  for il = 2 : nlayers
    tracersX(:,il) = [Xinit(1:end/2)+nx*dlayer(il);Xinit(end/2+1:end)+ny*dlayer(il)];
  end

  tracersXstore(:,:,idx) = tracersX;

  tracers.N = Nup;
  tracers.nv = nlayers-1;
  tracers.X = tracersX(:,2:nlayers);
 
  if 0
  figure(1); clf;
  plot(Xinit(1:end/2),Xinit(end/2+1:end),'k','linewidth',2)
  hold on
  plot(tracersX(1:end/2,:),tracersX(end/2+1:end,:),'ro','markersize',10)
  axis equal
  pause
  end

  % derivatives
  [Ben, Ten, Div] = vesicle.computeDerivs;

  % SLP 
  G = op.stokesSLmatrix(vesicle);

  % Build LHS and RHS
  LHS = (Div*G*Ten);
  RHS = -Div*(G*(Ben*Xinit));

  % solve for tension
  ten = LHS\RHS;

  % Calculate traction jump
  tracJump = vesicle.tracJump(Xinit,ten);

  % Velocity induced by traction jump
  selfVel = G*tracJump;
  kernel = @op.exactStokesSL;
  kernelDirect = @op.exactStokesSL;
  SLP = @(X) op.exactStokesSLdiag(vesicle,G,X);
  
  % Get the near zone
  [~,NearV2T] = vesicle.getZone(tracers,2);

  VelOnGrid = op.nearSingInt(vesicle,tracJump,SLP,[],NearV2T,kernel,kernelDirect,tracers,false,false);

  nearVelocity(:,:,idx) = [selfVel VelOnGrid];
  
  if 0
  figure(1);
  quiver(tracersX(1:end/2,:),tracersX(end/2+1:end,:),nearVelocity(1:end/2,:,idx),nearVelocity(end/2+1:end,:,idx))
  pause
  end

  idx = idx + 1;

  tend = toc(tstart);
  disp(['took ' num2str(tend) ' seconds'])
end

fileName = ['./output/nearFieldData/nearFieldData_' num2str(iset) '.mat']; 
save(fileName,'nSamplesInSet','XstandStore','nearVelocity','Nup','tracersXstore','-v7.3')
