clear; clc;
addpath ../src/
addpath ../examples/
oc = curve;
disp('Vesicle suspension in DLD')

% DLD parameters
%-------------------------------------------------------------------------
setBCmatFree = ~true; % true
prams.farField = 'rotDLD';
pillarCross = 'circlePillar'; % or trianglePillar
Dx = 10; % gap size in x-dir. (micron), if left [], determined based on wbox
Dy = 10; % gap size in y-dir. (micron)
prams.folderName = './output/polyDisp_VF20Cap75_circPost/';
wbox = 25; % square box size around a pillar (including Dx and Dy) (micron)
theta = atan(1/6);
prams.NbdExt = 1280; % # of points on exterior wall % 3072 for 6 rows
prams.NbdInt = 48; % # of points on posts % 64
prams.N = 32; % # of points per vesicle
prams.speed = 75; % Maximum velocity in the gap scale
prams.vesViscCont = 1; % viscosity contrast of the 1st type vesicle
prams.ves2ViscCont = 1; % viscosity contrast of the 2nd type vesicle
prams.ves2Concent = 0.1; % concentration of the 2nd type vesicle
prams.kappa = 1; % bending rigidity
prams.volFrac = 0.20; % volume fraction in a DLD lane
useFMM = false;

% FLAGS
%-------------------------------------------------------------------------
iplot = 0;
iJiggle = 1;
iuseNear = 1; % use wrong near-interactions (if =0, then neglect near-int)
iTrueNear = 0;
irepulsion = 0; % use repulsion, this is possible if iTrueNear == 1 or iuseNear == 1

% PARAMETERS, TOOLS
%-------------------------------------------------------------------------
exactSolveFreq = 0; % solve exactly at every [] time steps
prams.Th = 1.5; % time horizon
prams.nv = 1; prams.totnv = 1;
prams.fmm = useFMM; % use FMM for ves2ves
prams.fmmDLP = useFMM; % use FMM for ves2walls
prams.dt = 1E-4; % time step size

% NO NEED TO CHANGE ANYTHING BELOW THE LINE
% VESICLES and DLD DEVICE:
% -------------------------------------------------------------------------
scale = 0.05;
% 1) DLD Device
if strcmp(pillarCross,'circlePillar')
  prams.Dpost = 15*scale; % diameter of a post in micron * scale
  prams.Dpostx = prams.Dpost;
  prams.Dposty = prams.Dpost;
  prams.Dy    = Dy*scale; % vertical separation in micron * scale
  prams.Dx    = Dx*scale; % horizontal separation in micron * scale
  Xpost = [prams.Dpost/2*cos(-(0:prams.NbdInt-1)'/prams.NbdInt*2*pi);...
        prams.Dpost/2*sin(-(0:prams.NbdInt-1)'/prams.NbdInt*2*pi)];  
else
  load(pillarCross) % already scaled
  Xpost = [interpft(Xpost(1:end/2),prams.NbdInt);...
      interpft(Xpost(end/2+1:end),prams.NbdInt)];
  xup = interpft(Xpost(1:end/2),512); prams.Dpostx = max(xup)-min(xup);
  yup = interpft(Xpost(end/2+1:end),512); prams.Dposty = max(yup)-min(yup);
  if isempty(Dx)
    prams.Dx = scale*wbox-Dpostx;
    prams.Dy = scale*wbox-Dposty;
  else
    prams.Dx = scale*Dx;
    prams.Dy = scale*Dy;
  end
  prams.Dpost = [];
end
% Period of a device
prams.periodN = (prams.Dposty+prams.Dy)/(tan(theta)*(prams.Dpostx+prams.Dx)); 
% Row-shift fraction
prams.epsilon = 1/prams.periodN; 
% number of obstacles in x direction  (number of columns)
prams.ncol = 5;
% number of obstacles in y direction (number of rows)
prams.nrow    = 4; 

[XwallsInt,XwallsExt,prams.Lext,radiusRight,centx1stCol,centy1stCol,...
    prams.gPer] = initializeDLD(prams,Xpost,pillarCross,oc,setBCmatFree);
prams.xfreeze = centx1stCol+4*(prams.Dpostx+prams.Dx);


prams.nvbdInt = size(XwallsInt,2);
prams.nvbdExt = 1;
prams.xrange = min(XwallsExt(1:end/2));
% ------------------------
% -------------------------------------------------------------------------
% NET PARAMETERS
%-------------------------------------------------------------------------
Nnet = 256; % num. points
nComp = 16; % number of components for networks except relaxation problem
nCompRelax = 32; % number of PCA components for relaxation problem
% # of modes for M*vinf's network, inv(DivGT)*Div*vinf uses the same
nVelModes = 24; activeModes = [(1:nVelModes/2)';(Nnet-nVelModes/2+1:Nnet)'];
% number of modes for inv(DivGT)*(DivGB)x
nTenModes = 32; tenPredModes = [(1:nTenModes/2)';(Nnet-nTenModes/2+1:Nnet)'];
% load PCA basis vectors and column mean (evects, colMeans) for Nnet = 256
load necessaryMatFiles/pcaCoeffsBasis1step.mat
% load reference vesicle, then use it to randomly generate initial vesicles
load relaxedVesIC
dnn = dnnToolsDLD(X,XwallsInt,XwallsExt,prams);
tt = dnn.tt; dnn.oc = oc; dnn.useNear = iuseNear; dnn.useTrueNear = iTrueNear;
dnn.repulsion = irepulsion;
%-------------------------------------------------------------------------
disp(['Streaming vesicles in DLD, total #: ' num2str(prams.totnv) ', Kb = ' ...
    num2str(prams.kappa) ', Speed = ' num2str(prams.speed)])
%-------------------------------------------------------------------------
% LOAD NETWORKS (SEE THE END OF THE CODE FOR THESE FUNCTIONS)
% -------------------------------------------------------------------------
dnn.nComp = nComp;  dnn.nCompRelax = nCompRelax;
% LOAD PCA Network for Relaxation Problem
dnn = loadPCAnet4RelaxFiles(dnn);
% LOAD FFT Based Network for M*Vinf 
dnn = loadFCnet4AdvectFiles(nVelModes,activeModes,dnn);
% LOAD FFT Based Network for inv(DivGT)DivGBx
dnn = loadTenNetFiles(nTenModes,tenPredModes,dnn);
% LOAD FFT Based Network for inv(DivGT)*Div*vinf 
dnn = loadTenVnetFiles(nVelModes,activeModes,dnn);
% save PCA matrices 
dnn.colMeans = colMeans; dnn.evects = evects; 
% -------------------------------------------------------------------------

%--------------------------------------------------------------------------
% RANDOMLY INITIALIZE VESICLES
% -------------------------------------------------------------------------
[X,area0,len0,prams] = initializeVesicles(X,prams,XwallsInt,XwallsExt,tt,oc);
initCys = zeros(prams.nv,1);
for k = 1 : prams.nv
  initCys(k) = mean(X(end/2+1:end,k));
end
finalCys = zeros(prams.nv,1);

prams.initnv = prams.nv;

% Standardize vesicle RA = 0.9 and precompute the matrices 
[XstandEll,scalingEll,rotateEll,transEll,sortIdxEll] = ...
    dnn.standardizationStep(X(:,1),Nnet);
vesicle = capsules(XstandEll,[],[],1,1,prams.kappa);
vesicle.setUpRate([]);
opEll = poten(Nnet);
GEll = opEll.stokesSLmatrix(vesicle);
[BenEll,TenEll,DivEll] = vesicle.computeDerivs;
dnn.invTenMatEll = inv(DivEll*GEll*TenEll)*DivEll;
dnn.selfBendMatEll = inv(DivEll*GEll*TenEll)*(DivEll*GEll*(-BenEll))*prams.kappa;
dnn.Mell = GEll*TenEll*((DivEll*GEll*TenEll)\eye(Nnet))*DivEll;
dnn.relaxMatEll = (eye(2*Nnet)-prams.kappa*prams.dt*(-GEll*BenEll+dnn.Mell*GEll*BenEll))\eye(2*Nnet);

% load initConfigTest
% X = Xhist;
% [~,area0,len0] = oc.geomProp(X);
% prams.nv = size(X,2);
% prams.viscCont = ones(prams.nv,1);

Xfrozen = []; frozen_ids = [];
rem_ids = [1:prams.nv]'; nputLater = 0;
putBackLimit = ceil((prams.periodN+1)/4); 

disp(['Volume fraction ' num2str(prams.volFrac) ' leads to ' num2str(prams.nv) ' vesicles'])

%figure(1);clf;hold on;
%plot(XwallsExt(1:end/2),XwallsExt(end/2+1:end),'k','linewidth',2)
%plot(XwallsInt(1:end/2,:),XwallsInt(end/2+1:end,:),'k','linewidth',2)
%plot(X(1:end/2,:),X(end/2+1:end,:),'r','linewidth',2)
%axis equal
%pause

nPutBacks = 0; nZigZags = 0;
nves2 = ceil(prams.nv*prams.ves2Concent); % number of 2nd type of vesicles
ves2Ids = randperm(prams.nv,nves2); % global ids of 2nd type vesicles
% keep the viscosity contrasts of all cells to be streamed
allVesViscConts = ones(prams.nv,1)*prams.vesViscCont;
allVesViscConts(ves2Ids) = prams.ves2ViscCont;
% assign viscosity contrasts of initialized vesicles
prams.viscCont = allVesViscConts;   

% If using repulsion, find the scale
% -------------------------------------------------------------------------
if irepulsion
  vesicle = capsules(X(:,1),[],[],prams.kappa,1,true);
  vesicle.setUpRate();
  dnn.minDist = 0.4;
  dnn.repStrength = vesicle.repulStrengthScale(dnn.minDist,tt,prams.speed);
end    
% -------------------------------------------------------------------------

% INITIALLY TAKE SMALL TIME STEPS WITH IMPLICIT SOLVER
% ------------------------------------------------------------------------
% necessary to find correct initial tension, density, rotlet and stokeslet
tt.dt = 1E-5; sig = zeros(prams.N,prams.nv); 
etaInt = zeros(2*prams.NbdInt,prams.nvbdInt); 
etaExt = zeros(2*prams.NbdExt,1); RS = zeros(3,prams.nvbdInt+1);
for iter = 1 : 2
  vesicle = capsules(X,[],[],prams.kappa,ones(prams.nv,1),1); vesicle.setUpRate();
  [X,sig,etaInt,etaExt,RS] = tt.timeStepSimpleDiffDisc(X,sig,etaInt,etaExt,...
      RS,ones(prams.nv,1),dnn.wallsInt,dnn.wallsExt,vesicle);
end
tt.dt = prams.dt;
% ------------------------------------------------------------------------

% INITIALIZE MATRICES AND COUNTERS
% ------------------------------------------------------------------------  
time = (0:prams.dt:prams.Th)'; ntime = numel(time);
Xhist = X;
sigStore = sig;

fileName = [prams.folderName 'tStep1']; it = 1;
Dx = prams.Dx; Dy = prams.Dy; Dpostx = prams.Dpostx; Dposty = prams.Dposty;
speed = prams.speed; kappa = prams.kappa; volFrac = prams.volFrac;
save(fileName,'it','Xhist','Xfrozen','sigStore','XwallsInt','XwallsExt',...
    'rem_ids','frozen_ids','allVesViscConts','time','Dx','Dy','Dpostx',...
    'Dposty','kappa','speed','volFrac')
clear Dx Dy Dpostx Dposty speed kappa volFrac

ncountDNN = 0;
ncountExct = 0;
% ------------------------------------------------------------------------
% TIME STEPPING
for it = 2 : ntime
  disp('********************************************') 
  disp([num2str(it) 'th (/' num2str(ntime) ...
    ') time step, time: ' num2str(time(it))])

  
  if rem(it,exactSolveFreq) == 0       
    % EXACT SOLVE  
    disp('Taking an exact time step...'); tStart = tic;  
    ncountExct = ncountExct + 1;
    vesicle = capsules(Xhist,[],[],prams.kappa,ones(prams.nv,1),1);
    vesicle.setUpRate();
    [Xhist,sigStore,etaInt,etaExt,RS,...
        iter,iflag] = tt.timeStepSimpleDiffDisc(Xhist,...
        sigStore,etaInt,etaExt,RS,ones(prams.nv,1),...
        dnn.wallsInt,dnn.wallsExt,vesicle);
    
  else 
    % NETWORK SOLVE 
    ncountDNN = ncountDNN+1;
    disp('Taking a step with DNNs...');  tStart = tic;
    [Xhist,sigStore,etaInt,etaExt,RS] = dnn.DNNsolveMixed(...
        Xhist,sigStore,etaInt,etaExt,RS,Nnet);
   

    % check if there is any self-intersecting shape
    for k = 1 : prams.nv
      [xIntersect,~,~] = oc.selfintersect(Xhist(:,k));
      if ~isempty(xIntersect)
        disp('New vesicle shape is self-intersecting!!!')
      end
    end

    if ~isempty(xIntersect) 
      disp('Taking an exact time step...')
      ncountDNN = ncountDNN-1;
      ncountExct = ncountExct+1;
      vesicle = capsules(Xhist,[],[],prams.kappa,ones(prams.nv,1),1);
      vesicle.setUpRate();
      [Xhist,sigStore,etaInt,etaExt,RS,...
        iter,iflag] = tt.timeStepSimpleDiffDisc(Xhist,...
        sigStore,etaInt,etaExt,RS,ones(prams.nv,1),...
        dnn.wallsInt,dnn.wallsExt,vesicle);
    end % self-intersection occurs
    
  end % if rem(it,10) == 0

  % JIGGLING
  if iJiggle % jiggle vesicles pointwise if there is a near-collision
    Xhist = oc.fixCollisionsDLD(Xhist,XwallsInt);

    % Equally distribute points in arc-length b/c jiggling messes up points
    disp('Equally distributing points on arc-length...') 
    Xiter = Xhist;
    for iter = 1 : 5
      [Xiter,~,~] = oc.redistributeArcLength(Xiter);
    end

    % Fix misalignment in center and angle due to reparametrization
    Xhist = oc.alignCenterAngle(Xhist,Xiter);

  end
  
  % AREA-LENGTH CORRECTION
  disp('Correcting area-length errors...')
  [Xnew,ifail] = oc.correctAreaAndLength2(Xhist,area0,len0);
  
  if ifail
    disp('Error in AL cannot be corrected!!!')
  end
  Xhist = oc.alignCenterAngle(Xhist,Xnew);

  % check if shape is intersecting
  for k = 1 : prams.nv
    [xIntersect,~,~] = oc.selfintersect(Xhist(:,k));
    if ~isempty(xIntersect)
      disp('New vesicle shape is self-intersecting!!!')
    end
  end
  
  % Equally distribute points in arc-length
  disp('Equally distributing points on arc-length...')
  Xiter = Xhist;
  for iter = 1 : 5
    [Xiter,~,~] = oc.redistributeArcLength(Xiter);
  end

  % Fix misalignment in center and angle due to reparametrization
  Xhist = oc.alignCenterAngle(Xhist,Xiter);
  
  disp(['took ' num2str(toc(tStart)) ' seconds.'])
  
  % Compute error in area and length
  [~,area,len] = oc.geomProp(Xhist);
  errArea = max(abs(area-area0)./area0); errLen = max(abs(len-len0)./len0);
  errALPred = max(errArea,errLen);
  disp(['Error in area and length: ' num2str(errALPred)])   
  disp('********************************************') 
  disp(' ')
  
  
  % Then, put vesicles back and decide what to freeze 
  %[Xhist,sigStore,Xfrozen,rem_ids,frozen_ids,prams,nPutBacks,area0,len0] = putBackAndFreeze(...
  %      Xhist,sigStore,Xfrozen,prams,rem_ids,frozen_ids,allVesViscConts,nPutBacks,area0,len0);
    
  [Xhist,sigStore,Xfrozen,rem_ids,frozen_ids,prams,nPutBacks,...
    area0,len0,allVesViscConts,nputLater,nZigZags,initCys,finalCys] = freezeAndStream(oc,Xhist,sigStore,Xfrozen,prams,...
    rem_ids,frozen_ids,allVesViscConts,nPutBacks,area0,len0,nputLater,nZigZags,initCys,finalCys);

  prams.nv = numel(rem_ids);
  
  % First save current step, 
  fileName = [prams.folderName 'tStep' num2str(it)];
  save(fileName,'it','Xhist','Xfrozen','initCys','finalCys',...
      'sigStore','rem_ids','frozen_ids','nPutBacks','nputLater','nZigZags')
  if any(frozen_ids==1)
    disp('RA = 0.9 is frozen. COMPLETED!')
    break;
  end
  if numel(rem_ids) == 0 
    disp('all vesicles have zig-zagged. COMPLETED!')
    break;
  end
  if nPutBacks >= putBackLimit
    disp('Vesicle #1 has been put back at least # of period times. COMPLETED!')    
    break;
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xhist,sigStore,Xfrozen,rem_ids,frozen_ids,prams,nPutBacks,...
    area0,len0,allVesViscConts,nputLater,nZigZags,initCys,finalCys] = freezeAndStream(oc,X,sigma,Xfrozen,prams,...
    rem_ids,frozen_ids,allVesViscConts,nPutBacks,area0,len0,nputLater,nZigZags,initCys,finalCys)

% in this one we randomly put vesicle in the entrance if there is a vesicle 
% that reaches to the end of the device. If cannot put a new vesicle, then,
% freezes. In later time steps, it tries to put new vesicles.
    
% First check if the first vesicle (RA = 0.9) needs to be translated back
Xlarge = zeros(size(X));
for k = 1 : prams.nv
  cx = mean(X(1:end/2,k)); cy = mean(X(end/2+1:end,k));
  Xlarge(:,k) = [1.1*(X(1:end/2,k)-cx)+cx;1.1*(X(end/2+1:end,k)-cy)+cy];
end
vesicle = capsules(Xlarge,[],[],[],[],1);
vesicle.setUpRate();
NearV2V = vesicle.getZone([],1);
zoneV2V = NearV2V.zone;

pillarGrad = prams.epsilon*(prams.Dposty+prams.Dy)/(prams.Dpostx+prams.Dx);

% if vesicle is beyond xfreeze, then move it back to 1st gap
% HERE, ASSUMING THE INITIAL AREA WILL BE EMPTY FOR RA = 0.9

% ids of the vesicles to be frozen and to be remained
will_freeze=[]; will_remain=[];

centxn = mean(interpft(X(1:end/2,1),256));  
if centxn>=prams.xfreeze-prams.Dpostx/2
  centyn = mean(interpft(X(end/2+1:end,1),256));
  centxNew = centxn-4*(prams.Dx+prams.Dpostx);
  centyNew = centyn-pillarGrad*4*(prams.Dx+prams.Dpostx);
  X(:,1) = [X(1:end/2,1)-centxn+centxNew;X(end/2+1:end,1)-centyn+centyNew]; 
  nPutBacks = nPutBacks+1; 
    
  % check if there is any vesicle in near-zone of k, if so, move it back,
  % too
  K = (2:prams.nv);
  for k2 = K
    J = find(zoneV2V{k2}(:,1) == 1);
    if numel(J)~=0
      centyn = mean(interpft(X(end/2+1:end,k2),256));
      centxn = mean(interpft(X(1:end/2,k2),256));  
      centxNew = centxn-4*(prams.Dx+prams.Dpostx);
      centyNew = centyn-pillarGrad*4*(prams.Dx+prams.Dpostx);
      if min(X(1:end/2,k2)-centxn+centxNew)>prams.xrange
        X(:,k2) = [X(1:end/2,k2)-centxn+centxNew;X(end/2+1:end,k2)-centyn+centyNew];
        % if we move k2, then check its near zone and move vesicles close to
        % that too
        K2 = [(2:k2-1) (k2+1:prams.nv)];
        for k3 = K2
          J = find(zoneV2V{k3}(:,k2) == 1);
          if numel(J)~=0
            centyn = mean(interpft(X(end/2+1:end,k2),256));
            centxn = mean(interpft(X(1:end/2,k2),256));  
            centxNew = centxn-4*(prams.Dx+prams.Dpostx);
            centyNew = centyn-pillarGrad*4*(prams.Dx+prams.Dpostx);
            if min(X(1:end/2,k3)-centxn+centxNew)>prams.xrange
              X(:,k3) = [X(1:end/2,k3)-centxn+centxNew;X(end/2+1:end,k3)-centyn+centyNew];    
            else
              will_freeze = [will_freeze;k3];
              nputLater = nputLater + 1;
            end 
          end
        end % k3 = K2
      else % if we do not move k2, freeze it, also freeze those near to it
       will_freeze = [will_freeze;k2];
       nputLater = nputLater + 1;
       K2 = [(2:k2-1) (k2+1:prams.nv)];
       for k3 = K2
         J = find(zoneV2V{k3}(:,k2) == 1);
         if numel(J)~=0
           will_freeze = [will_freeze;k3];
           nputLater = nputLater + 1;
         end
       end % for k3 = K2
      end % if we move k2 near to RA = 0.9
    end % numel(J)~=0
  end % k2 = K
end

idOfLast = max([rem_ids;frozen_ids]); % id of last vesicle

% regions where a new vesicle can be introduced
xrange = prams.xrange + [0 prams.Dx+prams.Dpostx];
xrangePlace = prams.xrange + (prams.Dx+prams.Dpostx)/2;
yranges(1,:) = prams.Dy/2+prams.Dposty+[0.05*prams.Dy 0.95*prams.Dy]; % first lane
yranges(2,:) = -1.5*prams.Dy-prams.Dposty+[0.05*prams.Dy 0.95*prams.Dy]; % bottom lane
yranges(3,:) = [-0.45 0.45]*prams.Dy; % middle lane

% find vesicles' centers
centxS = mean(interpft(X(1:end/2,:),256),1)';
centyS = mean(interpft(X(end/2+1:end,:),256),1)';

% if we frozen some before, and need to seed again, randomly initialize vesicles
countPutLater = 0;
doItOneMore = true;
while (countPutLater<nputLater && doItOneMore)
  refVesId = randi([2,prams.nv]);
  Xref = X(:,refVesId);
  Xref = [Xref(1:end/2)-mean(Xref(1:end/2));Xref(end/2+1:end)-mean(Xref(end/2+1:end))];
  IA = oc.getIncAngle(Xref);
  xref = Xref(1:end/2)*cos(-IA)-Xref(end/2+1:end)*sin(-IA);
  yref = Xref(1:end/2)*sin(-IA)+Xref(end/2+1:end)*cos(-IA);
  Xref = [xref;yref];
  whichLane = randperm(2,2); % only in the bottom or top lanes
  vesPlaced = false; iter = 1;
  while (~vesPlaced && iter <= 2)
    % check if there is any vesicle in entrance of that lane
    vesExistsInlet = isempty(find(centyS >= yranges(whichLane(iter),1) & ...
          centyS <= yranges(whichLane(iter),2) & centxS >= xrange(1) & centxS <= xrange(2),1));
    if vesExistsInlet % then there is no vesicle in the inlet
      centxNew = xrangePlace;
      centyNew = (yranges(whichLane(iter),2)+yranges(whichLane(iter),1))*0.5;
      Xnew = [Xref(1:end/2)+centxNew; Xref(end/2+1:end)+centyNew];  
      vesPlaced = true; 
      countPutLater = countPutLater + 1;
      % update related information
      X(:,end+1) = Xnew;
      rem_ids = [rem_ids;idOfLast+1];
      area0(end+1) = area0(refVesId); len0(end+1) = len0(refVesId);
      prams.nv = prams.nv + 1; 
      allVesViscConts(end+1) = allVesViscConts(randi(numel(allVesViscConts)));
      prams.viscCont(end+1) = allVesViscConts(end);
      sigma(:,end+1) = sigma(:,refVesId);
      centxS(end+1) = centxNew; centyS(end+1) = centyNew;
      initCys(end+1) = centyNew;
      finalCys(end+1) = centyNew;
    else
      iter = iter + 1;
    end
  end % end while      
  doItOneMore = vesPlaced; % if vesicle is placed, then do it one more time, if we need to do so
end % nv <= prams.initnv
% update the number of vesicles to be put later
nputLater = nputLater - countPutLater;

% Get the near zone again
for k = 1 : prams.nv
  cx = mean(X(1:end/2,k)); cy = mean(X(end/2+1:end,k));
  Xlarge(:,k) = [1.1*(X(1:end/2,k)-cx)+cx;1.1*(X(end/2+1:end,k)-cy)+cy];
end
vesicle = capsules(Xlarge,[],[],[],[],1);
vesicle.setUpRate();
NearV2V = vesicle.getZone([],1);
zoneV2V = NearV2V.zone;


% if vesicle is beyond xfreeze, then move it back to
for k = 2 : prams.nv
  centxn = mean(interpft(X(1:end/2,k),256));  
  if centxn>=prams.xfreeze
    K = [(2:k-1) (k+1:prams.nv)];  
    centyn = mean(interpft(X(end/2+1:end,k),256));
    centxNew = centxn-4*(prams.Dx+prams.Dpostx);
    centyNew = centyn-pillarGrad*4*(prams.Dx+prams.Dpostx);
    vesExistsInlet = isempty(find(centyS(K) >= centyNew-0.5*prams.Dy & ...
        centyS(K) <= centyNew+0.5*prams.Dy & centxS(K) >= centxNew-prams.Dx/2-prams.Dpostx/2 & ...
        centxS(K) <= centxNew+prams.Dx/2+prams.Dpostx/2,1));
    % if there is no vesicle in the inlet  
    if vesExistsInlet && min(X(1:end/2,k)-centxn+centxNew)>prams.xrange
      X(:,k) = [X(1:end/2,k)-centxn+centxNew;X(end/2+1:end,k)-centyn+centyNew];
      % check if there is any vesicle in near-zone of k, if so, move it back,
      % too
      K = [(2:k-1) (k+1:prams.nv)];
      for k2 = K
        J = find(zoneV2V{k2}(:,k) == 1);
        if numel(J)~=0
          centyn = mean(interpft(X(end/2+1:end,k2),256));
          centxn = mean(interpft(X(1:end/2,k2),256));  
          centxNew = centxn-4*(prams.Dx+prams.Dpostx);
          centyNew = centyn-pillarGrad*4*(prams.Dx+prams.Dpostx);
          if min(X(1:end/2,k2)-centxn+centxNew)>prams.xrange
            X(:,k2) = [X(1:end/2,k2)-centxn+centxNew;X(end/2+1:end,k2)-centyn+centyNew];  
          else
            nputLater = nputLater + 1;
            will_freeze = [will_freeze; k2];
          end
        end
      end
    else % if there is a vesicle where the new one needs to be placed (as well as its near zone)
      nputLater = nputLater + 1; 
      will_freeze = [will_freeze; k]; 
      K = [(2:k-1) (k+1:prams.nv)];
      for k2 = K
        J = find(zoneV2V{k2}(:,k) == 1);
        if numel(J)~=0
          nputLater = nputLater + 1;
          will_freeze = [will_freeze; k2];
        end
      end
    end % if vesExistsInlet
  end % if centxn >= prams.xfreeze 
end % for k = 2 : prams.nv

for k = 1 : prams.nv
  finalCys(rem_ids(k)) = mean(X(end/2+1:end,k));
end

% Loop and check which ones zig-zagged and need to be frozen
yFreeze1stVes = -0.5*(prams.Dy+prams.Dposty);
if mean(X(end/2+1:end,1)) <= yFreeze1stVes
  % tag the ones we need to freeze
  will_freeze = [will_freeze;1];
  message = 'Vesicle 1 is frozen';
  disp(message)
end % if mean    

% use different threshold for others
delLat = (prams.Dy+prams.Dposty)*prams.epsilon;
yFreeze = -1.5*prams.Dy-2*prams.Dposty+4*delLat;
for iv = 2:prams.nv
  if mean(X(end/2+1:end,iv)) <= yFreeze
    % tag the ones we need to freeze
    will_freeze = [will_freeze;iv];
    nZigZags = nZigZags + 1;
    message = ['Vesicle ' num2str(rem_ids(iv),'%d') ' is frozen'];
    disp(message)
  end % if mean    
end % for iv

% save the frozen ones and remained ones
if~isempty(will_freeze)
  will_freeze = unique(will_freeze);
  for k = 1 : prams.nv
    if ~any(will_freeze == k)
      will_remain = [will_remain; k];
    end
  end
  % then there are vesicles to be frozen  
  Xfrozen = [Xfrozen X(:,will_freeze)];
  frozen_ids = [frozen_ids;rem_ids(will_freeze)]; 
else
  will_remain = [1:numel(X(1,:))]'; 
end

rem_ids = rem_ids(will_remain);
Xhist = X(:,will_remain);
sigStore = sigma(:,will_remain);
prams.viscCont = allVesViscConts(rem_ids);
prams.nv = numel(rem_ids);
area0 = area0(will_remain);
len0 = len0(will_remain);

end % freezeAndStream
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xhist,sigStore,Xfrozen,rem_ids,frozen_ids,prams,nPutBacks,area0,len0] = putBackAndFreeze(...
        X,sigma,Xfrozen,prams,rem_ids,frozen_ids,allVesViscConts,nPutBacks,area0,len0)

nv = size(X,2);
Xlarge = zeros(size(X));
for k = 1 : nv
  cx = mean(X(1:end/2,k)); cy = mean(X(end/2+1:end,k));
  Xlarge(:,k) = [1.5*(X(1:end/2,k)-cx)+cx;1.5*(X(end/2+1:end,k)-cy)+cy];
end
vesicle = capsules(Xlarge,[],[],[],[],1);
vesicle.setUpRate();
NearV2V = vesicle.getZone([],1);
zoneV2V = NearV2V.zone;

% if vesicle is beyond xfreeze, then move it back to
pillarGrad = prams.epsilon*(prams.Dposty+prams.Dy)/(prams.Dpostx+prams.Dx);
delLat = (prams.Dy+prams.Dposty)*prams.epsilon;
for k = 1 : nv
  centxn = mean(interpft(X(1:end/2,k),256));  
  if centxn>=prams.xfreeze
    centyn = mean(interpft(X(end/2+1:end,k),256));
    centxNew = centxn-4*(prams.Dx+prams.Dpostx);
    centyNew = centyn-pillarGrad*4*(prams.Dx+prams.Dpostx);
    X(:,k) = [X(1:end/2,k)-centxn+centxNew;X(end/2+1:end,k)-centyn+centyNew];
    if k == 1 % count nPutBacks for RA = 0.9 only
      nPutBacks = nPutBacks+1; 
    end
    
    % check if there is any vesicle in near-zone of k, if so, move it back,
    % too
    K = [(1:k-1) (k+1:nv)];
    for k2 = K
      J = find(zoneV2V{k2}(:,k) == 1);
      if numel(J)~=0
        centyn = mean(interpft(X(end/2+1:end,k2),256));
        centxn = mean(interpft(X(1:end/2,k2),256));  
        centxNew = centxn-4*(prams.Dx+prams.Dpostx);
        centyNew = centyn-pillarGrad*4*(prams.Dx+prams.Dpostx);
        X(:,k2) = [X(1:end/2,k2)-centxn+centxNew;X(end/2+1:end,k2)-centyn+centyNew];
        if k2 == 1
          nPutBacks = nPutBacks+1; 
        end
      end
    end
    
  end
end

yFreeze = -1.5*(prams.Dy+prams.Dposty)+4*delLat;

% ids of the vesicles to be frozen and to be remained
will_freeze=[]; will_remain=[];

% Loop over all the active vesicles
for iv = 1:size(X,2)

  if mean(X(end/2+1:end,iv)) <= yFreeze
    % tag the ones we need to freeze
    will_freeze = [will_freeze;iv];

    message = ['Vesicle ' num2str(rem_ids(iv),'%d') ' is frozen'];
    disp(message)
  else
    % tag the ones we need to keep
    will_remain = [will_remain;iv];
  end % if mean    
end % for iv

% save the frozen ones and remained ones
if~isempty(will_freeze)
  % then there are vesicles to be frozen  
  Xfrozen = [Xfrozen X(:,will_freeze)];
  frozen_ids = [frozen_ids;rem_ids(will_freeze)];  
end
rem_ids = rem_ids(will_remain);
Xhist = X(:,will_remain);
sigStore = sigma(:,will_remain);
prams.viscCont = allVesViscConts(rem_ids);
prams.nv = numel(rem_ids);
area0 = area0(will_remain);
len0 = len0(will_remain);


end % freezeAndStream

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,area0,len0,prams] = initializeVesicles(Xref,prams,XwallsInt,XwallsExt,tt,oc)

delLat = (prams.Dy+prams.Dposty)*prams.epsilon;

[~,area0,len0] = oc.geomProp(Xref);
x0 = interpft(Xref(1:end/2),prams.N);
y0 = interpft(Xref(end/2+1:end),prams.N);

% Left and right x-coordinate of the lanes (at the very left point 
% of the first pillar) + some buffer
%xrange(1) = min(XwallsInt(1:end/2,1))+1.5*prams.Dpostx+prams.Dx; % left point
xrange(1) = min(XwallsInt(1:end/2,1))+0.5*(prams.Dpostx+prams.Dx); % left point
xrange(2) = prams.xfreeze-0.5*(prams.Dpostx+prams.Dx); % right point

% initialize vesicles in the second half of the columns, then 
% put them back to the first half, so that they do not (?) collide when we
% put them back. Or find another way to impose periodicity
% maybe randomly put them (around their original position) but check if
% there is any collision at every time putting back
% nColsEmpty = ceil((ceil(prams.periodN)+2)/2);
% xrange(1) = min(XwallsInt(1:end/2,1)) + nColsEmpty*(prams.Dpostx+prams.Dx)-prams.Dx;
% xrange(2) = prams.xfreeze-prams.Dpostx/2;

% build objects for walls
wallsExt = capsules(XwallsExt,[],[],0,0,true);
wallsInt = capsules(XwallsInt,[],[],0,0,true);

% fill in 3 lanes (each is a parallelogram)
areaGeom = 3*(xrange(2)-xrange(1))*(0.95*prams.Dy);

% number of vesicles to fill 3 lanes
prams.nv = ceil(prams.volFrac*areaGeom/area0)+1; 

% Region info (1: above, 2: middle, 3: below)
regionyLB(1) = prams.Dy/2+prams.Dposty;
regionyLB(2) = -prams.Dy/2;
regionyLB(3) = -1.5*prams.Dy-prams.Dposty;
nvInRegs = zeros(3,1);

X = zeros(2*prams.N,prams.nv);

% place the RA = 0.9 particle
Xell = oc.ellipse(prams.N,0.9);
[~,~,lEll] = oc.geomProp(Xell);
Xell = Xell/lEll;

X(:,1) = [Xell(1:end/2)+min(XwallsInt(1:end/2,1))+1.5*prams.Dpostx+prams.Dx;...
  Xell(end/2+1:end)+delLat];

XLarge = X;

% counter
k = 2;

while k <= prams.nv

% randomly pick x-location for a center of vesicle
cx = xrange(1) + (xrange(2)-xrange(1))*rand;

% we need to pick cy now, but decide on the lane, so
% find the region which has the minimum number of vesicles 
[~,regID] = min(nvInRegs);

yBot = (cx-min(XwallsInt(1:end/2,1))-prams.Dpostx/2)*prams.epsilon+regionyLB(regID);
yTop = yBot+prams.Dy; 
yrange = [yBot+0.05*abs(yBot); yTop-0.05*abs(yTop)];

cy = yrange(1) + (yrange(2)-yrange(1))*rand;

phi = -pi/4+pi/2*rand;
% potential vesicle

xpot = cx + x0*cos(phi) + y0*sin(phi);
ypot = cy - x0*sin(phi) + y0*cos(phi);

xpotLarge = cx + 1.1*(x0*cos(phi) + y0*sin(phi));
ypotLarge = cy - 1.1*(x0*sin(phi) - y0*cos(phi));

accept = true; % tentatively accept the vesicle

% create capsule with the potential vesicle
vesicle = capsules([xpotLarge;ypotLarge],[],[],[],[],true);

[~,NearV2WInt] = vesicle.getZone(wallsInt,3);
[~,NearV2WExt] = vesicle.getZone(wallsExt,3);


[~,icollisionWallExt] = vesicle.collision(wallsExt,[],NearV2WExt,...
    tt.fmm,tt.op);
[~,icollisionWallInt] = vesicle.collision(wallsInt,[],NearV2WInt,...
    tt.fmm,tt.op);
if icollisionWallExt || icollisionWallInt
  accept = false;

  message = ['Vesicle crossed the outer wall.'];
  disp(message)
  % at least one of the vesicles's points is outside of one
  % of the solid wall components
end

if ~(cx <= xrange(2) && cx>= xrange(1) && cy <= yrange(2) ...
        && cy >= yrange(1))
    accept = false;

    message = ['Vesicle was outside the domain.'];
    disp(message)
end
% reject vesicle if it is outside the domain

if accept 
  % if vesicle is not outside of physical walls, accept it as
  % a potential new vesicle.  It will be kept as long as it intersects
  % no other vesicles  
  X(:,k) = [xpot;ypot];
  XLarge(:,k) = [xpotLarge;ypotLarge];

  % create an object with the current configuration of vesicles
  vesicle = capsules(XLarge(:,1:k),[],[],0,0,true);  


  % see if vesicles have crossed using collision detection code
  % Can be used with or without the fmm    
  NearV2V = vesicle.getZone([],1);
  icollisionVes = vesicle.collision(...
      [],NearV2V,[],tt.fmm,tt.op);

  if icollisionVes
    X(:,k) = 0;
    XLarge(:,k) = 0;
    message = ['Vesicles crossed.'];
    disp(message)
    % if they've crossed, reject it
  else
    nvInRegs(regID) = nvInRegs(regID) + 1;  
    k = k + 1;
    message = [num2str(prams.nv-k+1,'%d') ' vesicles left to fill the domain'];
    disp(message)
    % if they haven't crossed, increase the total number of vesicles
  end
end % if accept

end % while k <= 3*nv 

[~,area0,len0] = oc.geomProp(X);

end % initializeVesicles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [XwallsInt,XwallsExt,L,radiusRight,centx1stCol,centy1stCol,gPer] = ...
    initializeDLD(prams,Xpost,pillarCross,oc,setBCmatFree)
xup = interpft(Xpost(1:end/2),512); 
yup = interpft(Xpost(end/2+1:end),512);
% Move shape such that the rectangle around the shape is centered
x = Xpost(1:end/2); y = Xpost(end/2+1:end);
x2 = x-(max(xup)-prams.Dpostx/2);
y2 = y-(max(yup)-prams.Dposty/2);
Xpost = [x2;y2];

radiusLeft = prams.Dpostx/2; radiusRight = radiusLeft;

if ~strcmp(pillarCross,'circlePillar')
% Find where the smallest gap size is, place vesicles there and scale the
x2up = interpft(x2,512);
y2up = interpft(y2,512);

% look at the points in the horizontal direction and find the largest size
% across the shape in y-direction
maxGsize = -1E+4;
sort_x2 = sort(x2up);
for ij = 1 : numel(x2up)
  ids = find(abs(x2up-sort_x2(ij))<=1e-2);
  ys = y2up(ids);
  gapSize = max(ys)-min(ys);

  if gapSize>=maxGsize 
    maxGsize = gapSize;
    [~,idx] = max(y2up(ids));
    maxDist2Left = x2up(ids(idx))+prams.Dpostx/2;  
  end
end
else
maxDist2Left = prams.Dpost/2;    
end

% BUILD A PERIODIC DLD
H = prams.nrow*(prams.Dy+prams.Dposty);
L = prams.ncol*(prams.Dx+prams.Dpostx);

% Exterior long tube
a = L/2; b = H/2; order = 20;
t = (0:prams.NbdExt-1)'*2*pi/prams.NbdExt;
r = (cos(t).^order + sin(t).^order).^(-1/order);
x = a*r.*cos(t); y = b*r.*sin(t);

% ROTATE
pillarGrad = prams.epsilon*(prams.Dposty+prams.Dy)/(prams.Dpostx+prams.Dx);
y = y+pillarGrad*(x+L/2-prams.Dx/2-maxDist2Left);

XwallsExt = [x;y];

% Row shift  
delLat = (prams.Dy+prams.Dposty)*prams.epsilon;
% Place horizontally the obstacles in the first column  
centy1stCol = linspace(-H/2+prams.Dposty/2+prams.Dy/2,H/2-prams.Dposty/2-prams.Dy/2,prams.nrow)';
centx = linspace(-L/2+prams.Dpostx/2+prams.Dx/2,L/2-prams.Dpostx/2-prams.Dx/2,prams.ncol);

% center of the first column (x-coordinate)
centx1stCol = centx(1);

% Place vertically all the obstacles
centx = centx(ones(prams.nrow,1),:);

% Shift centers of the obstacles
delLatVect = [0:prams.ncol-1]*delLat;
centy = delLatVect(ones(prams.nrow,1),:) + centy1stCol(:,ones(1,prams.ncol));

% Place the obstacles
XwallsInt = zeros(2*prams.NbdInt,prams.nrow*prams.ncol);
for iwall = 1 : prams.nrow*prams.ncol
XwallsInt(:,iwall) = [Xpost(1:end/2)+centx(iwall);...
  Xpost(end/2+1:end)+centy(iwall)];
end


% BUILD A LARGE DLD, THEN COMPUTE BCs of a PERIODIC DLD
optionsWhole.farField = 'DLD';

optionsWhole.farFieldSpeed = prams.speed;
optionsWhole.confined = true;
optionsWhole.diffDiscWalls = true; % walls are discretized by different Nbd
optionsWhole.antiAlias = true;
pramsWhole.runName = 'WholeDLD';
optionsWhole.logFile  = './output/wholeDLD.log';
% Name of binary data file for storing vesicle information
optionsWhole.dataFile = './output/wholeDLDData.bin';
pramsWhole.periodN = prams.periodN;
pramsWhole.ncol = 9;
pramsWhole.nrow = 8; % 10
pramsWhole.epsilon = 1/prams.periodN;
pramsWhole.Dpost = prams.Dpost;
pramsWhole.Dx = prams.Dx;
pramsWhole.Dy = prams.Dy;
pramsWhole.NbdInt = prams.NbdInt;
pramsWhole.NbdExt = 3328; % 3584
pramsWhole.gmresTol = 1e-7;  % tolerance for gmres

optionsWhole.fastDirect = false;
optionsWhole.fmm = false;  % fmm for single-layer potentials
optionsWhole.fmmDLP = false; % fmm for double-layer potentials
optionsWhole.matFreeWalls = false; % W2W interactions are done without a matrix
optionsWhole.HODLRforW2W  = false;  % Use HODLR to compress wall2wall interactions
[optionsWhole,pramsWhole] = initVes2D(optionsWhole,pramsWhole);

[XwallsInt_whole,XwallsExt_whole,Lwhole,Hwhole,~,~,~,...
    pramsWhole.Dpostx,pramsWhole.Dposty] = ...
    oc.initConfigDLD('wholeDLD',pramsWhole.NbdInt,pramsWhole.NbdExt,...
    Xpost,pramsWhole.Dpost,pramsWhole.epsilon,pramsWhole.periodN,...
    pramsWhole.Dx,pramsWhole.Dy,pramsWhole.nrow,pramsWhole.ncol);
pramsWhole.nvbdInt = size(XwallsInt_whole,2);
pramsWhole.nvbdExt = 1;
pramsWhole.nvbd = pramsWhole.nvbdInt + 1;

% center of the first post in the periodic device
xcentper = centx1stCol; 
ycentper = -1.5*(prams.Dy+prams.Dposty);

% this center is supposed to move to
xcentwhole = -Lwhole/2+3*prams.Dx+3.5*prams.Dpostx-maxDist2Left;
ycentwhole = -2.5*(prams.Dposty+prams.Dy) + ...
    pillarGrad*2*(prams.Dpostx+prams.Dx);

% the amount we need to move is 
dispx = xcentwhole-xcentper;
dispy = ycentwhole-ycentper;

% target points (exterior wall of periodic model)
Xtra = [XwallsExt(1:end/2)+dispx;XwallsExt(end/2+1:end)+dispy];

% first compute velocity at the middle and scale the maximum velocity
gapX = -Lwhole/2+[2;3;4;5;6;7;8]*(pramsWhole.Dpostx+pramsWhole.Dx);
gapY = [1;2;3;4;5;6;7]*pillarGrad*(pramsWhole.Dpostx+pramsWhole.Dx);
gapY1 = gapY-2*prams.Dposty-2*prams.Dy;
gapY2 = gapY-prams.Dposty-prams.Dy;

  
% DEBUG
idebug = 0;
if idebug
  figure(1); clf; hold on;
  plot([XwallsExt_whole(1:end/2);XwallsExt_whole(1)],[XwallsExt_whole(end/2+1:end);XwallsExt_whole(end/2+1)],'k','linewidth',2)
  plot([XwallsInt_whole(1:end/2,:);XwallsInt_whole(1,:)],[XwallsInt_whole(end/2+1:end,:);XwallsInt_whole(end/2+1,:)],'k','linewidth',2)
  plot([Xtra(1:end/2);Xtra(1)],[Xtra(end/2+1:end);Xtra(end/2+1)],'r-o','markersize',1,'markerfacecolor','r','linewidth',2)
  plot([XwallsInt(1:end/2,:);XwallsInt(1,:)]+dispx,[XwallsInt(end/2+1:end,:);XwallsInt(end/2+1,:)]+dispy,'r','linewidth',2)
  plot(gapX,gapY,'bo','markersize',6,'markerfacecolor','b')
  axis equal
  pause
end

% construct tstep and compute velocity
omWhole = monitor([],optionsWhole,pramsWhole);
ttWhole = tstep(optionsWhole,pramsWhole,omWhole);

if ~setBCmatFree
  [vel,~,~] = ttWhole.computeVelFieldNoVes(pramsWhole,[],XwallsInt_whole,...
      XwallsExt_whole,[gapX;gapX;gapX;gapY1;gapY;gapY2],0);
else
  [vel,~,~] = ttWhole.computeVelFieldMatFree(pramsWhole,XwallsInt_whole,...
      XwallsExt_whole,[gapX;gapX;gapX;gapY1;gapY;gapY2],[],0);  
end
% Scale the maximum velocity input to match the max. velocity at the gap
scalingVel = optionsWhole.farFieldSpeed/mean(vel(1:end/2));
optionsWhole.farFieldSpeed = optionsWhole.farFieldSpeed*scalingVel;

% Set the BC for the whole device again
ttWhole.farField = @(X,Xint) ttWhole.bgFlow(X,optionsWhole.farField,...
    'Speed',optionsWhole.farFieldSpeed,...
    'intWalls',XwallsInt_whole,'nrow',pramsWhole.nrow,...
    'ncol',pramsWhole.ncol,'Dpostx',pramsWhole.Dpostx,...
    'Dposty',pramsWhole.Dposty,'GapX',pramsWhole.Dx,'GapY',pramsWhole.Dy,...
    'epsilon',pramsWhole.epsilon);
% Now compute the velocity on the periodic device's exterior wall
% this is fast
if ~setBCmatFree
  [gPer,~,~] = ttWhole.computeVelFieldNoVes(pramsWhole,[],XwallsInt_whole,...
      XwallsExt_whole,Xtra,0);
else
  [gPer,~,~] = ttWhole.computeVelFieldMatFree(pramsWhole,XwallsInt_whole,...
      XwallsExt_whole,Xtra,scalingVel,0);  
end

% Remove large matrices of the whole device from memory
clear ttWhole;
end % initializeDLD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dnn = loadPCAnet4RelaxFiles(dnn)
% network for solving relaxation problem for coupled tension and X

if dnn.kappa == 1 % KAPPA = 1 NETWORKS (Dts = 1e-5, 5e-5, 1e-4, 5e-4, 1e-3)
if dnn.dt == 1e-5
  file = './networks/fcPCArelaxN256Dt1E5nModes';
elseif dnn.dt == 5e-5
  file = './networks/fcPCArelaxN256Dt5E5nModes';    
elseif dnn.dt == 1e-4
  file = './networks/fcPCArelaxN256nModes';
elseif dnn.dt == 5e-4
  file = './networks/fcPCArelaxN256Dt5E4nModes';
elseif dnn.dt == 1e-3
  file = './networks/fcPCArelaxN256Dt1E3nModes';
else
  disp('there is no network trained for the chosen time step size, stop!')
  pause
end

% normalized output
load([file '1to16_fcXlarge_tstepFCNet_w1step.mat'])
% load('./networks/fcPCArelax_fcXlarge_tstepFCNet_w1step.mat')
dnn.bendNets{1} = net; 
dnn.muChan_bend(1) = muChan1; 
dnn.sdevChan_bend(1) = sdevChan1; 
dnn.scale_bend(1) = scale; 
dnn.offset_bend(1) = offset;
if dnn.nCompRelax > 16
load([file '17to32_fcXlarge_tstepFCNet_w1step.mat'])
dnn.bendNets{2} = net; 
dnn.muChan_bend(2) = muChan1; 
dnn.sdevChan_bend(2) = sdevChan1; 
dnn.scale_bend(2) = scale; 
dnn.offset_bend(2) = offset;
end
if dnn.nCompRelax > 32
load([file '33to48_fcXlarge_tstepFCNet_w1step.mat'])
dnn.bendNets{3} = net; 
dnn.muChan_bend(3) = muChan1; 
dnn.sdevChan_bend(3) = sdevChan1; 
dnn.scale_bend(3) = scale; 
dnn.offset_bend(3) = offset;
end
if dnn.nCompRelax > 48
load([file '49to64_fcXlarge_tstepFCNet_w1step.mat'])
dnn.bendNets{4} = net; 
dnn.muChan_bend(4) = muChan1; 
dnn.sdevChan_bend(4) = sdevChan1; 
dnn.scale_bend(4) = scale; 
dnn.offset_bend(4) = offset;
end

elseif dnn.kappa == 1e-1 % KAPPA = 1E-1 NETWORKS (Dt = 1E-4 and 5E-5 ONLY)

if dnn.dt == 1e-4 
file = './networks/fcPCArelaxN256Dt1E4Kb1E1nModes';
elseif dnn.dt == 5e-5
file = './networks/fcPCArelaxN256Dt5E5Kb1E1nModes';
end

% normalized output
load([file '1to16_fcXlarge_tstepFCNet_w1step.mat'])
dnn.bendNets{1} = net; 
dnn.muChan_bend(1) = muChan1; 
dnn.sdevChan_bend(1) = sdevChan1; 
dnn.scale_bend(1) = scale; 
dnn.offset_bend(1) = offset;
if dnn.nCompRelax > 16
load([file '17to32_fcXlarge_tstepFCNet_w1step.mat'])
dnn.bendNets{2} = net; 
dnn.muChan_bend(2) = muChan1; 
dnn.sdevChan_bend(2) = sdevChan1; 
dnn.scale_bend(2) = scale; 
dnn.offset_bend(2) = offset;
end
if dnn.nCompRelax > 32
load([file '33to48_fcXlarge_tstepFCNet_w1step.mat'])
dnn.bendNets{3} = net; 
dnn.muChan_bend(3) = muChan1; 
dnn.sdevChan_bend(3) = sdevChan1; 
dnn.scale_bend(3) = scale; 
dnn.offset_bend(3) = offset;
end
if dnn.nCompRelax > 48
load([file '49to64_fcXlarge_tstepFCNet_w1step.mat'])
dnn.bendNets{4} = net; 
dnn.muChan_bend(4) = muChan1; 
dnn.sdevChan_bend(4) = sdevChan1; 
dnn.scale_bend(4) = scale; 
dnn.offset_bend(4) = offset;
end    
    
    
elseif dnn.kappa == 1e-2 % KAPPA = 1E-2 NETWORKS (Dt = 1E-4 ONLY)

file = './networks/fcPCArelaxN256Dt1E4Kb1E2nModes';

% normalized output
load([file '1to16_fcXlarge_tstepFCNet_w1step.mat'])
dnn.bendNets{1} = net; 
dnn.muChan_bend(1) = muChan1; 
dnn.sdevChan_bend(1) = sdevChan1; 
dnn.scale_bend(1) = scale; 
dnn.offset_bend(1) = offset;
if dnn.nCompRelax > 16
load([file '17to32_fcXlarge_tstepFCNet_w1step.mat'])
dnn.bendNets{2} = net; 
dnn.muChan_bend(2) = muChan1; 
dnn.sdevChan_bend(2) = sdevChan1; 
dnn.scale_bend(2) = scale; 
dnn.offset_bend(2) = offset;
end
if dnn.nCompRelax > 32
load([file '33to48_fcXlarge_tstepFCNet_w1step.mat'])
dnn.bendNets{3} = net; 
dnn.muChan_bend(3) = muChan1; 
dnn.sdevChan_bend(3) = sdevChan1; 
dnn.scale_bend(3) = scale; 
dnn.offset_bend(3) = offset;
end
if dnn.nCompRelax > 48
load([file '49to64_fcXlarge_tstepFCNet_w1step.mat'])
dnn.bendNets{4} = net; 
dnn.muChan_bend(4) = muChan1; 
dnn.sdevChan_bend(4) = sdevChan1; 
dnn.scale_bend(4) = scale; 
dnn.offset_bend(4) = offset;
end    
    
    
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dnn = loadFCnet4AdvectFiles(nmodes,activeModes,dnn)
% network for M acting on vback

% initialize mu,sdev. for PCA net, we do not output any of them.
muChan1 = []; sdevChan1 = []; scale = []; offset = []; outputSize = [];
%netFilesFolder = './networks/fcNetVelocityPredFilesPCAin/velPredPCAin_mode';
netFilesFolder = './networks/n256Mtimes24modesFFTNets/velPredPCAin_mode';

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
load('./networks/fcPCAtenMaTonBendN256_32modes_FCNet_w1step.mat')
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
netFilesFolder = './networks/n256tenMatTimes24modesFFTNets/velPredPCAin_mode';

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
