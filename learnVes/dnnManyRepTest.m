function dnnManyRepTest(repStrength, repScale)
addpath ../src/
addpath ../examples/
disp('Multiple vesicles with background fluid, NNA scheme')

% FLAGS
%-------------------------------------------------------------------------
variableKbDt = true;
interpOrder = 3;
prams.bgFlow = 'couette'; % 'rotation' or 'couette' (w/ solid boundaries)
prams.speed = 100; % 70 for rotation, 100 for couette 
iDNNlike = 0; % whether solve exactly DNNlike or use DNNs
iplot = 0;
iJiggle = 1;
iuseNear = 0; % use wrong near-interactions (if =0, then neglect near-int)
iTrueNear = 1;
irepulsion = 1; % use repulsion, this is possible if iTrueNear == 1
saveFreq = 1;
downTo = 16;
% PARAMETERS, TOOLS
%-------------------------------------------------------------------------
errTol = 1e-2;
maxDt = 1e-4;
exactSolveFreq = 0; % solve exactly at every [] time steps
prams.Th = 1.5; % time horizon
prams.N = 32; % num. points for true solve in DNN scheme
prams.Nexact = 128;
prams.nv = 70; %(24 for VF = 0.1, 47 for VF = 0.2) num. of vesicles
prams.fmm = ~false; % use FMM for ves2ves
prams.fmmDLP = ~false; % use FMM for ves2walls
prams.kappa = 1;
prams.interpOrder = interpOrder;
prams.dtExact = 1E-4;
prams.dt = 1E-4; % time step size
if strcmp(prams.bgFlow,'couette')
  prams.Nbd = 1536;
  prams.nvbd = 2;
else
  prams.Nbd = 0;
  prams.nvbd = 0;
end
oc = curve;
Th = prams.Th; N = prams.N; nv = prams.nv; dt = prams.dt; 
Nbd = prams.Nbd; nvbd = prams.nvbd; bgFlow = prams.bgFlow; speed = prams.speed;

% net parameters
Nnet = 256; % num. points
nComp = 16; % number of components for networks except relaxation problem
nCompRelax = 32; % number of PCA components for relaxation problem
% # of modes for M*vinf's network, inv(DivGT)*Div*vinf uses the same
nVelModes = 24; activeModes = [(1:nVelModes/2)';(Nnet-nVelModes/2+1:Nnet)'];
% number of modes for inv(DivGT)*(DivGB)x
nTenModes = 32; tenPredModes = [(1:nTenModes/2)';(Nnet-nTenModes/2+1:Nnet)'];
% load PCA basis vectors and column mean (evects, colMeans) for Nnet = 256
% load necessaryMatFiles/pcaCoeffsBasis1step.mat
load necessaryMatFiles/pcaBasisNewest.mat
%-------------------------------------------------------------------------
disp(['Flow: ' prams.bgFlow ', N = ' num2str(N) ', nv = ' num2str(nv) ...
    ', Th = ' num2str(Th)])
%-------------------------------------------------------------------------

% VESICLES and WALLS:
% -------------------------------------------------------------------------
initType = 3; % 1: VesID, 2: vesShape, 3: initialVesFile, 4: randomly fill
vesID = 88201;  % 88201 from the library
vesShape = 'ellipse'; % 'ellipse' or 'curly' 
initialVesFile = 'VF30_70VesIC';%'VF30_70VesIC'; % from another run % probVF20IC: fails in 2it
%initialVesFile = 'deRefIssueIC';
%initialVesFile = 'failConfig';
cent = []; % centers of vesicles, if empty, randomly assigned 
thet = []; % angular positions of vesicles, if empty, randomly assigned
IA = []; % initial inclination angles of vesicles
irandInit = 1; % perturb center, angle and IA of vesicles?
volFrac = 0.3;
[X,Xwalls,area0,len0] = initializeVesiclesAndWalls(vesID,...
    vesShape,initialVesFile,cent,thet,IA,initType,bgFlow,N,nv,Nbd,nvbd,irandInit,oc,[],[]);
% -------------------------------------------------------------------------

if iDNNlike
  solveType = 'DNNlike';
else
  solveType = 'DNN';
end

fileName = ['./output/repulsion/repStr' num2str(repStrength) '_repScale' num2str(repScale) 'VF30.bin'];
fid = fopen(fileName,'w');
output = [N;nv;Nbd;nvbd];
fwrite(fid,output,'double');
x = X(1:end/2,:); y = X(end/2+1:end,:); xwalls = Xwalls(1:end/2,:); ywalls = Xwalls(end/2+1:end,:);
output = [x(:);y(:);xwalls(:);ywalls(:)];
fwrite(fid,output,'double');
fclose(fid);

% BUILD DNN CLASS
% -------------------------------------------------------------------------
dnn = dnnToolsManyVesNoLoop(X,Xwalls,prams);
tt = dnn.tt; dnn.oc = oc; dnn.useNear = iuseNear; dnn.useTrueNear = iTrueNear;
ttExact = dnn.ttExact;
dnn.repulsion = irepulsion; dnn.variableKbDt = variableKbDt;
if irepulsion
  vesicle = capsules(X(:,1),[],[],prams.kappa,1,true);
  vesicle.setUpRate();
  dnn.minDist = repScale; % dist was 1 for 1E6
  [dnn.repStrength, dnn.repLenScale] = vesicle.repulStrengthScaleSimple(dnn.minDist,tt,prams.speed);
  dnn.repStrength = repStrength;
  disp(['Repulsion strength = '  num2str(dnn.repStrength)])
end    
% -------------------------------------------------------------------------

% UNCOMMENT TO GENERATE AN IC FOR DENSE COUETTE
%[X,Xwalls,area0,len0] = initializeVesiclesAndWalls(vesID,...
%    vesShape,initialVesFile,cent,thet,IA,4,bgFlow,N,nv,Nbd,nvbd,irandInit,oc,volFrac,tt);
%figure(1);clf;hold on;
%plot(Xwalls(1:end/2,:),Xwalls(end/2+1:end,:),'k');
%plot(X(1:end/2,:),X(end/2+1:end,:),'r')
%axis equal
%pause


% LOAD NETWORKS (SEE THE END OF THE CODE FOR THESE FUNCTIONS)
% -------------------------------------------------------------------------
if ~iDNNlike
dnn.nComp = nComp;  dnn.nCompRelax = nCompRelax;
% LOAD PCA Network for Relaxation Problem
dnn = loadAllPCAnets4Relax(dnn);
% LOAD FFT Based Network for M*Vinf 
dnn = loadFCnet4AdvectFiles(nVelModes,activeModes,dnn);
% LOAD FFT Based Network for inv(DivGT)DivGBx
dnn = loadTenNetFiles(nTenModes,tenPredModes,dnn);
% LOAD FFT Based Network for inv(DivGT)*Div*vinf 
dnn = loadTenVnetFiles(nVelModes,activeModes,dnn);
% save PCA matrices 
dnn.colMeans = colMeans; dnn.evects = evects; 
end
% -------------------------------------------------------------------------

% INITIALLY TAKE SMALL TIME STEPS WITH IMPLICIT SOLVER
% ------------------------------------------------------------------------
% necessary to find correct initial tension, density, rotlet and stokeslet
% sig = zeros(prams.N,nv); eta = zeros(2*Nbd, nvbd); RS = zeros(3,nvbd);
sig = zeros(prams.Nexact,nv); eta = zeros(2*Nbd,nvbd); RS = zeros(3,nvbd);
if numel(X(:,1))/2 ~= prams.Nexact
  X = [interpft(X(1:end/2,:),prams.Nexact); interpft(X(end/2+1:end,:),prams.Nexact)];
end
for iter = 1 : 2
 vesicle = capsules(X,[],[],prams.kappa,ones(nv,1),1); vesicle.setUpRate();
 [X,sig,eta,RS] = ttExact.timeStepSimple(X,sig,eta,RS,ones(nv,1),dnn.walls,vesicle);
end
if numel(X(:,1))/2 ~= prams.N
  X = [interpft(X(1:end/2,:),prams.N); interpft(X(end/2+1:end,:),prams.N)];
  sig = interpft(sig,prams.N);
end
% ------------------------------------------------------------------------

% INITIALIZE MATRICES AND COUNTERS
% ------------------------------------------------------------------------
time = [0];
iDNNsolved = [1];
Xhist = X; sigStore = sig; etaStore = zeros(2*Nbd,nvbd); RSstore = zeros(3,nvbd);
errALPred = zeros(1,1);
if Nbd ~= 0
  etaStore(:,:) = eta; RSstore(:,:) = RS;
end

ncountCNN = 0;
ncountExct = 0;
% ------------------------------------------------------------------------

writeData(fileName,Xhist,sigStore,etaStore,RSstore,time(end),ncountCNN,ncountExct,0);

% TIME STEPPING
it = 1;
while time(end) < prams.Th
  disp('********************************************') 
  disp([num2str(it) 'th time step, time: ' num2str(time(it))])

  shouldSolveExact = 0;
  justSolvedExact = 0;
  nRelax = 0;
  if rem(it,exactSolveFreq) == 0       
    % EXACT SOLVE  
    disp('Taking an exact time step...'); tStart = tic;  
    ncountExct = ncountExct + 1;
    iDNNsolved(it) = 0;
    X = Xhist; sig = sigStore;
    if numel(Xhist(:,1))/2 ~= prams.Nexact 
      X = [interpft(Xhist(1:end/2,:),prams.Nexact); interpft(Xhist(end/2+1:end,:),prams.Nexact)];
      sig = interpft(sigStore,prams.Nexact);
    end
    [Xreparam,niter] = oc.reparametrize(X,[],6,50);
    X = oc.alignCenterAngle(X,Xreparam);
    vesicle = capsules(X,[],[],prams.kappa,ones(nv,1),1);
    vesicle.setUpRate();
    [Xnew,signew,etanew,RSnew,iter,iflag] = ttExact.timeStepSimple(X,sig,...
        etaStore,RSstore,ones(nv,1),dnn.walls,vesicle);
    justSolvedExact = 1;
    if numel(Xhist(:,1))/2 ~= prams.Nexact 
      Xnew = [interpft(Xnew(1:end/2,:),prams.N); interpft(Xnew(end/2+1:end,:),prams.N)];
      signew = interpft(signew,prams.N);
    end
    
  else 
    % NETWORK SOLVE 
    ncountCNN = ncountCNN+1;
    iDNNsolved(it) = 0;
    disp('Taking a step with DNNs...');  tStart = tic;
    [Xnew,signew,etanew,RSnew,vback] = dnn.DNNsolve2(Xhist,sigStore,etaStore,RSstore,Nnet);
    % check if there is any self-intersecting shape
    for k = 1 : nv
      [xIntersect,~,~] = oc.selfintersect(Xnew(:,k));
      if ~isempty(xIntersect)
        disp('MLARM gives a shape that is self-intersecting, SOO')
        shouldSolveExact = 0;
        Xinter = Xhist(:,k);
        save intersectingVesShape Xinter
        % Relax that shape with its background velocity
        disp('Relaxing it for a few time step')
        nRelax = nRelax + 1;
        Xup = [interpft(Xhist(1:end/2,k),prams.Nexact); interpft(Xhist(end/2+1:end,k),prams.Nexact)];
        vup = [interpft(vback(1:end/2,k),prams.Nexact); interpft(vback(end/2+1:end,k),prams.Nexact)];
        for jts = 1 : 5
          vesicle = capsules(Xup, [], [], prams.kappa, ones(1,1), 1);
          vesicle.setUpRate();
          Xup = dnn.relaxExactSolve(vesicle, vup, dt/5, Xup, ttExact.op);
        end
        Xdown = [interpft(Xup(1:end/2),prams.N);interpft(Xup(end/2+1:end),prams.N)];
        [xIntersect,~,~] = oc.selfintersect(Xdown);
        if isempty(xIntersect)
          disp('Self-intersection is resolved') % this is exact solve for only the intersecting vesicle
          Xnew(:,k) = Xdown;
        else
          disp('Self-intersection cannot be resolved')
          save afterRelaxVesShape Xdown
          break;
        end
      end
      % OR SMOOTH OUT THAT SHAPE AND SOLVE WITH DNN AGAIN,
      % OR RELAX THAT SHAPE FOR A FEW SMALL TIME STEPS AND SOLVE WITH DNN
      % AGAIN
    end
    if shouldSolveExact 
      disp('TAKING AN EXACT TIME STEP...')
      ncountCNN = ncountCNN-1;
      ncountExct = ncountExct+1;
      iDNNsolved(it) = 0;
      X = Xhist; sig = sigStore;
      if numel(Xhist(:,1))/2 ~= prams.Nexact 
        X = [interpft(Xhist(1:end/2,:),prams.Nexact); interpft(Xhist(end/2+1:end,:),prams.Nexact)];
        sig = interpft(sigStore,prams.Nexact);
      end
      [Xreparam,niter] = oc.reparametrize(X,[],6,50);
      X = oc.alignCenterAngle(X,Xreparam);
      vesicle = capsules(X,[],[],prams.kappa,ones(nv,1),1);
      vesicle.setUpRate();
      [Xnew,signew,etanew,RSnew,iter,iflag] = ttExact.timeStepSimple(X,sig,etaStore,RSstore,ones(nv,1),dnn.walls,vesicle);
      if numel(Xhist(:,1))/2 ~= prams.Nexact 
        Xnew = [interpft(Xnew(1:end/2,:),prams.N); interpft(Xnew(end/2+1:end,:),prams.N)];
        signew = interpft(signew,prams.N);
      end
      justSolvedExact = 1;
    end % self-intersection occurs
    
  end % if rem(it,10) == 0
  

  % JIGGLING
  if iJiggle % jiggle vesicles pointwise if there is a near-collision
    Xnew = oc.fixCollisionsZinc(Xnew,Xwalls);
  end
  
  % AREA-LENGTH CORRECTION
  disp('Correcting area-length errors...')
  [Xnew2,ifail] = oc.correctAreaAndLength3(Xnew,area0,len0,downTo);
  
  % COLLISION CHECK
  Xup = [interpft(Xnew2(1:end/2,:),prams.Nexact);interpft(Xnew2(end/2+1:end,:),prams.Nexact)];
  vesicleProv = capsules(Xup, [], [], [], [], 1);
  vesicleProv.setUpRate(ttExact.op);
  collUpRate = vesicleProv.uprate;
  Nup = size(Xup,1)/2 * collUpRate;
  Xup =  [interpft(Xup(1:end/2,:),Nup);interpft(Xup(end/2+1:end,:),Nup)];
  vesicleUp = capsules(Xup,[],[],[],[],1);
  [NearV2V,NearV2W] = vesicleUp.getZone(dnn.walls,3);
  [icollisionVes,icollisionWall] = ...
    vesicleUp.collision(dnn.walls,NearV2V,NearV2W,ttExact.fmm,ttExact.op);
  
  if icollisionVes || icollisionWall
    disp('There is a collision')
    save collisionCase Xnew
    if ~justSolvedExact
      disp('Solving again, now with the implicit solver')
      ncountExct = ncountExct+1;
      ncountCNN = ncountCNN - 1;
      iDNNsolved(it) = 0;
      X = Xhist; sig = sigStore;
      if numel(Xhist(:,1))/2 ~= prams.Nexact 
        X = [interpft(Xhist(1:end/2,:),prams.Nexact); interpft(Xhist(end/2+1:end,:),prams.Nexact)];
        sig = interpft(sigStore,prams.Nexact);
      end
      [Xreparam,niter] = oc.reparametrize(X,[],6,50);
      X = oc.alignCenterAngle(X,Xreparam);
      vesicle = capsules(X,[],[],prams.kappa,ones(nv,1),1);
      vesicle.setUpRate();
      [Xnew,signew,etanew,RSnew,iter,iflag] = ttExact.timeStepSimple(X,sig,etaStore,RSstore,ones(nv,1),dnn.walls,vesicle);
      if numel(Xhist(:,1))/2 ~= prams.Nexact 
        Xnew = [interpft(Xnew(1:end/2,:),prams.N); interpft(Xnew(end/2+1:end,:),prams.N)];
        signew = interpft(signew,prams.N);
      end
      disp('Correcting area-length errors...')
      [Xnew2,ifail] = oc.correctAreaAndLength3(Xnew,area0,len0,downTo);
    end
  end
  Xnew = oc.alignCenterAngle(Xnew,Xnew2);
  
  Xup = [interpft(Xnew(1:end/2,:),prams.Nexact);interpft(Xnew(end/2+1:end,:),prams.Nexact)];
  vesicleProv = capsules(Xup, [], [], [], [], 1);
  vesicleProv.setUpRate(ttExact.op);
  collUpRate = vesicleProv.uprate;
  Nup = size(Xup,1)/2 * collUpRate;
  Xup =  [interpft(Xup(1:end/2,:),Nup);interpft(Xup(end/2+1:end,:),Nup)];
  vesicleUp = capsules(Xup,[],[],[],[],1);
  [NearV2V,NearV2W] = vesicleUp.getZone(dnn.walls,3);
  [icollisionVes,icollisionWall] = ...
    vesicleUp.collision(dnn.walls,NearV2V,NearV2W,ttExact.fmm,ttExact.op);

   
  disp(['took ' num2str(toc(tStart)) ' seconds.'])
  % Compute error in area and length
  [~,area,len] = oc.geomProp(Xnew);
  errArea = max(abs(area-area0)./area0); errLen = max(abs(len-len0)./len0);
  dtnew = prams.dt;
  
  
  if max(errArea, errLen) >= errTol || icollisionVes || icollisionWall
    dtnew = 0.5 * prams.dt; % Reject
    if dtnew > maxDt; dtnew = maxDt; end;
  elseif max(errArea, errLen) < 0.98 * errTol
    dtnew = 1.2 * prams.dt; % Accept
    if dtnew > maxDt; dtnew = maxDt; end;
    it = it + 1;
    Xhist = Xnew;
    sigStore = signew;
    etaStore = etanew;
    RSstore = RSnew;

    errALPred(it) = max(errArea,errLen);
    time(it) = time(it-1) + dtnew;  
  else % still accept but do not change the time step size
    it = it + 1;
    Xhist = Xnew;
    sigStore = signew;
    etaStore = etanew;
    RSstore = RSnew;
    errALPred(it) = max(errArea,errLen);
    time(it) = time(it-1) + dtnew;  
    if dtnew > maxDt; dtnew = maxDt; end;
  end
  disp(['Error in area and length: ' num2str(max(errArea, errLen))])   
  disp(['New time step size: ' num2str(dtnew)])
  disp('********************************************') 
  disp(' ')
  
  if rem(it,saveFreq) == 0
    disp(['Number of times a vesicle is relaxed isolated = ' num2str(nRelax)])  
    writeData(fileName,Xhist,sigStore,etaStore,RSstore,time(end),ncountCNN,ncountExct,iDNNsolved(end));
  end

  if dtnew ~= prams.dt
  tt.dt = dtnew;
  ttExact.dt = dtnew;
  prams.dt = dtnew;
  dnn.dt = dtnew;
  dnn.tt = tt;
  dnn.ttExact = ttExact;
  dnn = loadAllPCAnets4Relax(dnn);
  end
  
  if dtnew < 1e-6
    disp('Time step size is too small')
    break;
  end

  if iplot
  figure(1);clf;
  plot(Xwalls(1:end/2,:),Xwalls(end/2+1:end,:),'k')
  hold on;
  plot(Xhist(1:end/2,:),Xhist(end/2+1:end,:),'r','linewidth',2)
  axis equal
  pause(0.1)
  end
end

% Save data to a mat-file:
writeData(fileName,Xhist,sigStore,etaStore,RSstore,time(end),ncountCNN,ncountExct,iDNNsolved(end));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Xwalls,area0,len0] = initializeVesiclesAndWalls(vesID,...
    vesShape,initialVesFile,cent,thet,IA,initType,bgFlow,N,nv,Nbd,nvbd,irandInit,oc,volFrac,tt)
Xwalls = zeros(2*Nbd,nvbd);
if strcmp(bgFlow,'couette')
  thet = (0:Nbd-1)'*2*pi/Nbd;
  Xwalls = [ [2.2*cos(thet); 2.2*sin(thet)] [cos(-thet); sin(-thet)] ];
end

if initType == 1
  load ./necessaryMatFiles/errInPred100K_FCnetPCA.mat
  errAllPCA = errPCA;
  load /workspace/gokberk/X100KinitShapes.mat  
  X0 = Xstore(:,vesID);
  X0 = [interpft(X0(1:end/2),N);interpft(X0(end/2+1:end),N)];
  %errFFT = ceil(errAllFFT(vesID)*100);
  errPCA = ceil(errAllPCA(vesID)*100);
  disp(['VesID: ' num2str(vesID) ', Error in prediction: (PCA) ' ...
      num2str(errPCA) '% '])
elseif initType == 2 
  X0 = oc.initConfig(N,vesShape);
  [~,~,len] = oc.geomProp(X0);
  X0 = X0./len;
  disp(['VesShape: ' vesShape ])
elseif initType == 3
  disp('loading an IC file')
  load(initialVesFile)   
  if numel(X(:,1))/2 ~= N
    X = [interpft(X(1:end/2,:),N);interpft(X(end/2+1:end,:),N)];
  end 
  [~,area0,len0] = oc.geomProp(X);
elseif initType == 4
  X0 = oc.initConfig(N,'ellipse');
  [~,~,len] = oc.geomProp(X0);
  X0 = X0./len;  
end

if isempty(IA)
  if irandInit
    IA = 2*pi*rand(nv,1);
  else
    IA = zeros(nv,1);
  end
end

if initType < 3
  [~,area0,len0] = oc.geomProp(X0);
  vesRadius = sqrt(area0(1)/pi); % 0.1291
  area0 = ones(1,nv)*area0; len0 = ones(1,nv)*len0;

  X = zeros(2*N,nv);
  % Move center for rotation flow
  if strcmp(bgFlow,'rotation')
    if isempty(thet)
      thet = (0:pi/nv:pi*(nv-1)/nv)'; 
      if irandInit
        thet = 0.95*thet + 0.1*thet*rand(nv,1);
      end 
    end
    if isempty(cent)
      cent = 10*vesRadius;
    end
    
    for k = 1 : nv
      x = X0(1:end/2)*cos(IA(k))-X0(end/2+1:end)*sin(IA(k));
      y = X0(1:end/2)*sin(IA(k))+X0(end/2+1:end)*cos(IA(k));
      X(:,k) = [x+cent*cos(thet(k));y+cent*sin(thet(k))];
    end  
  end
  if strcmp(bgFlow,'couette')
    if isempty(thet)
      thet = (0:pi/nv:pi*(nv-1)/nv)'; 
      if irandInit
        thet = (thet-0.05) + 0.1*rand(nv,1);
      end
    end
    if isempty(cent)   
      if irandInit
        cent = 1.6 + 0.15*rand(nv,1);
      else
        cent = ones(nv,1)*1.6;
      end
    end
    for k = 1 : nv
      x = X0(1:end/2)*cos(IA(k))-X0(end/2+1:end)*sin(IA(k));
      y = X0(1:end/2)*sin(IA(k))+X0(end/2+1:end)*cos(IA(k));  
      X(:,k) = [x+cent(k)*cos(thet(k)); y+cent(k)*sin(thet(k))];   
    end
  end      
elseif initType == 4 % initType == 4 so randomly fill domain
   X = oc.fillCouetteArea(X0,Xwalls,volFrac,tt);
   [~,area0,len0] = oc.geomProp(X);
end % if initType



end % initializeVesiclesAndWalls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dnn = loadAllPCAnets4Relax(dnn)
interpOrder = dnn.interpOrder;
% list of networks for different Kb*Dt values
dnn.KbDts = 2.^(0:12)' * 1E-6;
files = [];
% fileName = './flow_networks/fcPCArelaxN256Dt';
fileName = './networks2/fcPCArelaxN256Dt';
for ifl = 1 : numel(dnn.KbDts)
  files{ifl} = [fileName num2str(dnn.KbDts(ifl)) 'Kb1nModes'];
end

% current flow's Kb*Dt
flowKbDt = dnn.kappa * dnn.dt;

[~, idx] = min(abs(dnn.KbDts-flowKbDt));
    
if idx == 1 || idx == numel(dnn.KbDts)
  disp('Extrapolation needed for the given bending stiffness and Dt, stop!')
  whichNets = [idx];
else
  % Choose 5 networks and we do 5th order Lagrange Interpolation 
  if idx == numel(dnn.KbDts)-1
    whichNets = (idx-(interpOrder-2):idx+1);
  elseif idx == 2
    whichNets = (idx-1:idx+(interpOrder-2));
  else
    whichNets = (idx-(interpOrder-1)/2:idx+(interpOrder-1)/2);
  end
%   whichNets = (1:numel(dnn.KbDts))';
end

% Load the networks needed for Lagrange interpolation
for k = 1 : numel(whichNets)
  if 1 %dnn.KbDts(k) >= 5E-4
  load([files{whichNets(k)} '1to16_fcXlarge_tstepFCNet_flow.mat'])
  else
  load([files{whichNets(k)} '1to16_fcMedium_tstepFCNet_flow.mat'])    
  end
  dnn.bendNets{k,1} = net; 
  dnn.muChan_bend(k,1) = muChan1; 
  dnn.sdevChan_bend(k,1) = sdevChan1; 
  dnn.scale_bend(k,1) = scale; 
  dnn.offset_bend(k,1) = offset;
  if dnn.nCompRelax > 16
  if 1 %dnn.KbDts(k) >= 5E-4
  load([files{whichNets(k)} '17to32_fcXlarge_tstepFCNet_flow.mat'])
  else
  load([files{whichNets(k)} '17to32_fcMedium_tstepFCNet_flow.mat']) 
  end
  dnn.bendNets{k,2} = net; 
  dnn.muChan_bend(k,2) = muChan1; 
  dnn.sdevChan_bend(k,2) = sdevChan1; 
  dnn.scale_bend(k,2) = scale; 
  dnn.offset_bend(k,2) = offset;
  end  
end
dnn.KbDts = dnn.KbDts(whichNets);

end % loadAllPCAnets4Relax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dnn = loadFCnet4AdvectFiles(nmodes,activeModes,dnn)
% network for M acting on vback

% initialize mu,sdev. for PCA net, we do not output any of them.
muChan1 = []; sdevChan1 = []; scale = []; offset = []; outputSize = [];
%netFilesFolder = './networks/fcNetVelocityPredFilesPCAin/velPredPCAin_mode';
netFilesFolder = './networks/NEWn256Mtimes24modesFFTNets/velPredPCAin_mode';

for imode = 1 : nmodes
  pmode = activeModes(imode);
  load([netFilesFolder num2str(pmode) '_net_w1step.mat'])
  FCnets{imode} = net;
end

dnn.MVnets = FCnets; dnn.muChan_MV = muChan1; dnn.sdevChan_MV = sdevChan1;
dnn.scale_MV = scale; dnn.offset_MV = offset; dnn.MVoutSize = outputSize;
dnn.nVelModes = nmodes; dnn.velActiveModes = activeModes;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dnn = loadTenNetFiles(nTenModes,tenPredModes,dnn)
% network for inverse of tension matrix on self bending

%load('./networks/fcPCAtenMatOnBendN256_FCNet_w1step.mat')
load('./networks/NEWfcPCAtenMaTonBendN256_32modes_FCNet_w1step.mat')
dnn.tenBendNets = net; dnn.muChan_tenBend = muChan1; 
dnn.sdevChan_tenBend = sdevChan1; dnn.scale_tenBend = scale; 
dnn.offset_tenBend = offset; dnn.nTenModes = nTenModes; 
dnn.tenPredModes = tenPredModes;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dnn = loadTenVnetFiles(nmodes,activeModes,dnn)
% network for inverse of tension matrix on vback

% initialize mu,sdev. for PCA net, we do not output any of them.
muChan1 = []; sdevChan1 = []; scale = []; offset = []; outputSize = [];
%netFilesFolder = './networks/n256tenMatTimesFFTNets/velPredPCAin_mode';
netFilesFolder = './networks/NEWn256tenMatTimes24modesFFTNets/velPredPCAin_mode';

for imode = 1 : nmodes
  pmode = activeModes(imode);
  load([netFilesFolder num2str(pmode) '_net_w1step.mat'])
  FCnets{imode} = net;
end

dnn.tenVnets = FCnets; dnn.muChan_tenV = muChan1; dnn.sdevChan_tenV = sdevChan1;
dnn.scale_tenV = scale; dnn.offset_tenV = offset; dnn.tenVoutSize = outputSize;
dnn.nTenVmodes = nmodes; dnn.TenVactiveModes = activeModes;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeData(filename,X,sigma,eta,RS,time,ncountNN,ncountExact,iDNNsolved)
x = X(1:end/2,:);
y = X(end/2+1:end,:);
output = [time;ncountNN;iDNNsolved;x(:);y(:);sigma(:)];
if ~isempty(eta)
etax = eta(1:end/2,:); etay = eta(end/2+1:end,:);
output = [output;etax(:);etay(:);RS(:)];    
end

fid = fopen(filename,'a');
fwrite(fid,output,'double');
fclose(fid);


end

