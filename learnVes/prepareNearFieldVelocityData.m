function prepareNearFieldVelocityData(iset,npar)

% Load Xstore containing nves number of vesicles 
% Xstore holds 2*N entries in its columns as [x;y] coordinates of nves
% number of vesicles
load /work2/03353/gokberk/frontera/X100KinitShapes.mat
oc = curve;

% Information about loaded data
nves = size(Xstore,2); % number of vesicles in the data set
Nstore = size(Xstore,1)/2; % number of discretization points on vesicles


% Upsampling maybe necessary
Nup = 256;


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
nearVelocity = zeros(2*Nup,4,nSamplesInSet);
tracersXstore = zeros(2*Nup,4,nSamplesInSet);
XstandStore = zeros(2*Nup,nSamplesInSet);

idx = 1;
for k = sum(nSamples(1:iset-1))+1 : sum(nSamples(1:iset))
  disp(['Vesicle #' num2str(idx) ' out of ' num2str(nSamples(k)) ' being processed...'])
  tstart = tic;

  % change the resolution
  Xinit = [interpft(Xstore(1:end/2,k),Nup); interpft(Xstore(end/2+1:end,k),Nup)]; 
  % Standardize if not done yet
  if ~istandard
    [Xinit,scaling,rotate,trans,sortIdx] = dnn.standardizationStep(Xinit,oc);
  end
  [~,area,len] = oc.geomProp(Xinit);
  
  % store the shape
  XstandStore(:,idx) = Xinit;
  
  % Generate the grid for velocity
  [~,tang] = oc.diffProp(Xinit);
  % get x and y components of normal vector at each point
  nx = tang(Nup+1:2*Nup);
  ny = -tang(1:Nup);

  % Points where velocity is calculated involve the points on vesicle
  tracersX = zeros(Nup, 4);
  tracersX(:,1) = Xinit;

  % initialize vesicle
  vesicle = capsules(Xinit, [], [], 1, 1, 1);
  vesicle.setUpRate();

  % Generate tracers
  h = vesicle.length/vesicle.N;  % arc-length spacing
  tracersX(:,2) = [Xinit(1:end/2)+nx*h/2; Xinit(end/2+1:end)+ny*h/2]; % at h/2
  tracersX(:,3) = [Xinit(1:end/2)+nx*h; Xinit(end/2+1:end)+ny*h]; % at h
  tracersX(:,4) = [Xinit(1:end/2)+nx*1.5*h; Xinit(end/2+1:end)+ny*1.5*h]; % at 3/2 * h
  
  tracersXstore(:,:,idx) = tracersX;

  tracers.N = numel(Xtra)/2;
  tracers.nv = 1;
  tracers.X = tracersX;
 

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
  idx = idx + 1;

  tend = toc(tstart);
  disp(['took ' num2str(tend) ' seconds'])
end

fileName = ['/work2/03353/gokberk/frontera/nearFieldData/Data_' num2str(iset) '.mat']; 
save(fileName,'nSamplesInSet','XstandStore','nearVelocity','Nup','tracersXstore','-v7.3')
