function computeStressTensor(fileName,runName,kappa,VC,farField,Ufar,skip,fmm)
% ONLY FOR CONFINED FLOWS (EITHER DLD OR ANY OTHER)
% Target points is the points on the vesicle
addpath ../
addpath ../../src/

% Load Data
[posx,posy,wallx,wally,extWallx,extWally,intWallx,intWally,...
time,N,nv] = loadFile(fileName);

% Check if there are walls with different discretizations
iDiffDiscWalls = 0; %i.e. if DLD flow or the other confined flows
options.diffDiscWalls = 0;
if isempty(wallx)
  iDiffDiscWalls = 1;
  options.diffDiscWalls = 1;
end

% Parameters and Options
prams.kappa = kappa; % bending stiffness
prams.viscCont = VC; % viscosity contrast
prams.N = N; % number of points per vesicle

% walls' parameters
prams.NbdInt = size(intWallx,1);
prams.NbdExt = size(extWallx,1);
prams.Nbd = size(wallx,1);
prams.nvbdInt = size(intWallx,2);
prams.nvbdExt = size(extWallx,2);
prams.nvbd = prams.nvbdInt+prams.nvbdExt+size(wallx,2);
prams.nv = size(posx,2);

prams.nrow = 12;
prams.ncol = 9;
prams.Dpostx = 1.5;
prams.Dposty = 1.5;
prams.Dx = 1;
prams.Dy = 1;
prams.epsilon = 1/6;

prams.gmresTol = 1e-8;
options.fmm = fmm;
options.fmmDLP = fmm;
options.matFreeWalls = true;

options.farField  = farField; % background velocity
options.farFieldSpeed = Ufar; % scaling of farField
options.confined  = true;        % confined
options.order     = 1;           % time stepping order
options.inextens  = 'method1';   % Inextensibility condition
options.vesves    = 'implicit';  
options.near      = true;        % near-singular integration
options.orderGL   = 2;    % ~ Gauss-Lobatto
options.nsdc      = 0;    % # of SDC corrections
options.antiAlias = true; % Anti-aliasing

prams.runName = runName;
options.logFile = ['./' prams.runName '.log'];
options.dataFile = ['./' prams.runName '_Data.bin']; 
[options,prams] = initVes2D(options,prams);

% Build structures
X = [posx(:,1);posy(:,1)];
om = monitor(X,options,prams);
tt = tstep(options,prams,om);
Xwalls = [wallx;wally];
XwallsInt = [intWallx;intWally];
XwallsExt = [extWallx;extWally];

% Build the walls and the preconditioner, too
if ~iDiffDiscWalls
  walls = tt.initialConfined(prams,Xwalls,[],[]);
else
  [~,walls_interior,walls_exterior] = tt.initialConfined(prams,[],...
      XwallsInt,XwallsExt);
end

% poten classes 
potWall    = tt.opWall;
potWallInt = tt.opWallInt;
potWallExt = tt.opWallExt;
potVes     = tt.op;

% If we skip some time steps:
timeSkipped = time(1:skip:end);
ntimeSkipped = numel(timeSkipped);

% storage matrices
sigmaStore = zeros(prams.N,ntimeSkipped); % for tension
interfacialForceStore = zeros(2*prams.N,ntimeSkipped);

% for stress(xx,xy) due to vesicle
stress1SelfStore = zeros(2*prams.N,ntimeSkipped); 
% for stress(xy,yy) due to vesicle
stress2SelfStore = stress1SelfStore;

% for stress(xx,xy) due to walls (or interior walls)
stress1Store = zeros(2*prams.N,ntimeSkipped);
% for stress(xy,yy) due to walls (or interior walls)
stress2Store = stress1Store;

% for stress(xx,xy) due to exterior wall  
stress1ExtStore = stress1Store;
% for stress(xy,yy) due to exterior wall
stress2ExtStore = stress2Store;    


totStress1 = stress1Store;
totStress2 = stress2Store;

% Compute tension, rotlets and Stokeslets, density
idx = 1;
for k = 1 :skip: numel(time)
  tTime = tic;
  
  % Current vesicle configuration
  X = [posx(:,:,k);posy(:,:,k)];
  
  % If you ask pressure around a vesicle, put Xnew scaled up as below
%   Xnew = 1.2*[X(1:end/2)-mean(X(1:end/2));X(end/2+1:end)-mean(X(end/2+1:end))];
%   Xnew = [Xnew(1:end/2)+mean(X(1:end/2));Xnew(end/2+1:end)+mean(X(end/2+1:end))];
  
  % build vesicle class and class for target points at which stress
  % computed
  vesicle = capsules(X,[],[],prams.kappa,prams.viscCont,1);
  vesicle.setUpRate(potVes);

  %stressTar = capsules(X,[],[],[],[],1);
  %stressTar.setUpRate(potVes);
  
  % Get the near structre
%   [~,NearV2T] = vesicle.getZone(stressTar,2);
  
%  if iDiffDiscWalls
%    [~,NearWint2T] = walls_interior.getZone(stressTar,2);
%    [~,NearWext2T] = walls_exterior.getZone(stressTar,2); 
%  else
%    [~,NearW2T] = walls.getZone(stressTar,2);
%  end

  % Compute tension, density and RS
  disp('Computing tension,density and RS')

  if ~iDiffDiscWalls
    [sigma,eta,~,~,RS,~,iter] = vesicle.computeSigAndEta(tt,walls,[],[]);
  else
    [sigma,~,etaInt,etaExt,RS,~,iter] = vesicle.computeSigAndEta(tt,[],...
        walls_interior,walls_exterior);
  end
  
%   vesicle.sig = sigma;
  f = vesicle.tracJump(X,sigma);
  
  interfacialForceStore(:,idx) = f;
  sigmaStore(:,idx) = sigma;
   
  % Compute stress
  if 0
  disp('Computing stress')
  
  if iDiffDiscWalls
    % Due to the pillars
    [stress1,stress2] = ...
      walls_interior.stressTensor(etaInt,RS(:,2:end),...
        stressTar,options.fmm,'DLP',potWallInt,NearWint2T);
    stress1Store(:,idx) = stress1;
    stress2Store(:,idx) = stress2;
    % Due to the exterior wall
    [stress1,stress2] = ...
      walls_exterior.stressTensor(etaExt,RS(:,1),...
        stressTar,options.fmm,'DLP',potWallExt,NearWext2T);
    stress1ExtStore(:,idx) = stress1;
    stress2ExtStore(:,idx) = stress2;
  else
    % Due to walls  
    [stress1,stress2] = walls.stressTensor(eta,RS,stressTar,options.fmm,...
        'DLP',potWall,NearW2T);
    stress1Store(:,idx) = stress1;
    stress2Store(:,idx) = stress2;
    
  end
  
  % Due to the vesicle itself
  [stress1,stress2] = ...
    vesicle.stressSelfTensor(f,'SLP',potVes);
  stress1SelfStore(:,idx) = stress1;
  stress2SelfStore(:,idx) = stress2;
  % if there is a viscosity contrast add that, too??
%   if VC~=1
%     [stress1,stress2] = ...
%       vesicle.stressSelfTensor(f,'DLP',potVes);
%     stress1SelfStore(:,idx) = stress1SelfStore(:,idx) + stress1;
%     stress2SelfStore(:,idx) = stress2SelfStore(:,idx) + stress2; 
%   end
  
  % Compute the total stress
  totStress1(:,idx) = stress1Store(:,idx) + stress1ExtStore(:,idx) + ...
      stress1SelfStore(:,idx);
  totStress2(:,idx) = stress2Store(:,idx) + stress2ExtStore(:,idx) + ...
      stress2SelfStore(:,idx);
  end

  idx = idx + 1;
  
  tTime = toc(tTime);
  message = [num2str(k) 'th time step out of ' num2str(numel(time)) ...
    ' time steps takes ' num2str(tTime,'%2.2d') ' seconds'];
  om.writeMessage(message,'%s\n')
end

fileExp = [runName '_Stress.mat'];
save(fileExp,'posx','posy','time','interfacialForceStore','skip','stress1Store','stress2Store',...
    'stress1ExtStore','stress2ExtStore','stress1SelfStore','stress2SelfStore',...
    'totStress1','totStress2','sigmaStore','wallx','wally','intWallx','intWally',...
    'extWallx','extWally','kappa','VC','iDiffDiscWalls');
end


