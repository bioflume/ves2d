classdef dnnClass

properties
oc
KbDts   
tt
dt
kappa
vinf
variableKbDt

confined
walls
NearW2W
wallsInt
wallsExt
NearWint2Wint
NearWint2Wext
NearWext2Wint

bendNets
muChan_bend
sdevChan_bend
scale_bend
offset_bend
MVnets
muChan_MV
sdevChan_MV
scale_MV
offset_MV
MVoutSize
tenBendNets
muChan_tenBend
sdevChan_tenBend
scale_tenBend
offset_tenBend
nTenModes
tenPredModes
tenVnets
muChan_tenV
sdevChan_tenV
scale_tenV
offset_tenV
tenVoutSize
nTenVmodes
TenVactiveModes
nVelModes
velActiveModes
nComp
nCompRelax
colMeans
evects
Nnet

nCircular
invTenMatCirc
selfBendMatCirc
Mcirc
relaxMatCirc

end

methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = dnnClass(X,Xwalls,XwallsInt,XwallsExt,prams,dnnOpts)
if nargin > 0
% Flow parameters    
o.confined = dnnOpts.confined;
o.kappa = prams.kappa;
o.dt = prams.dt;    
o.nCircular = prams.nCircular;
o.tt = o.buildTstep(X,prams);  
if o.confined
  o = o.initWalls(Xwalls,XwallsInt,XwallsExt);
  o.vinf = [];
else
  o.vinf = o.setBgFlow(prams.bgFlow,prams.speed);  
end % if confined

% Network options
o.Nnet = dnnOpts.Nnet;
o.nVelModes = dnnOpts.nVelModes;
o.nTenModes = dnnOpts.nTenModes;
o.nComp = dnnOpts.nComp;
o.nCompRelax = dnnOpts.nCompRelax;

% load networks
o = o.loadPCAnets4Relax;
o = o.loadFCnet4AdvectFiles;
o = o.loadTenNetFiles;
o = o.loadTenVnetFiles;

end % if nargin > 0
end % dnnClass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = loadPCAnets4Relax(o)

% list of networks for different Kb*Dt values
o.KbDts = [5e-7;1e-6;5e-6;1e-5;5e-5;1e-4;2.5E-4;5e-4;7.5E-4;1e-3;2.5E-3];

files{11} = './networks/fcPCArelaxN256Dt2p5E34Kb1nModes';
files{10} = './networks/fcPCArelaxN256Dt1E3Kb1nModes';
files{9} = './networks/fcPCArelaxN256Dt7p5E4Kb1nModes';
files{8} = './networks/fcPCArelaxN256Dt5E4Kb1nModes';
files{7} = './networks/fcPCArelaxN256Dt2p5E4Kb1nModes';
files{6} = './networks/fcPCArelaxN256Dt1E4Kb1nModes';
files{5} = './networks/fcPCArelaxN256Dt5E4Kb1E1nModes';    
files{4} = './networks/fcPCArelaxN256Dt1E4Kb1E1nModes';
files{3} = './networks/fcPCArelaxN256Dt5E5Kb1E1nModes';
files{2} = './networks/fcPCArelaxN256Dt1E4Kb1E2nModes';
files{1} = './networks/fcPCArelaxN256Dt5E5Kb1E2nModes';

% current flow's Kb*Dt
flowKbDt = o.kappa * o.dt;

[~, idx] = min(abs(o.KbDts-flowKbDt));
    
if idx == 1 || idx == numel(o.KbDts)
  disp('Extrapolation needed for the given bending stiffness and Dt, stop!')
  pause
else
  % Choose 5 networks and we do 5th order Lagrange Interpolation 
  if idx == numel(o.KbDts)-1
    whichNets = (idx-3:idx+1);
  elseif idx == 2
    whichNets = (idx-1:idx+3);
  else
    whichNets = (idx-2:idx+2);
  end
end

% if flowKbDt is equal to one of the trained values, then use that network
o.variableKbDt = true;
if any(o.KbDts == flowKbDt)
  whichNets = find(o.KbDts == flowKbDt);
  o.variableKbDt = false;
end

% Load the networks needed for Lagrange interpolation
for k = 1 : numel(whichNets)
  load([files{whichNets(k)} '1to16_fcXlarge_tstepFCNet_w1step.mat'])
  o.bendNets{k,1} = net; 
  o.muChan_bend(k,1) = muChan1; 
  o.sdevChan_bend(k,1) = sdevChan1; 
  o.scale_bend(k,1) = scale; 
  o.offset_bend(k,1) = offset;
  if o.nCompRelax > 16
  load([files{whichNets(k)} '17to32_fcXlarge_tstepFCNet_w1step.mat'])
  o.bendNets{k,2} = net; 
  o.muChan_bend(k,2) = muChan1; 
  o.sdevChan_bend(k,2) = sdevChan1; 
  o.scale_bend(k,2) = scale; 
  o.offset_bend(k,2) = offset;
  end  
end
o.KbDts = o.KbDts(whichNets);


end % loadPCAnets4Relax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = loadFCnet4AdvectFiles(o)

% Fourier modes
o.velActiveModes = [(1:o.nVelModes/2)';(o.Nnet-o.nVelModes/2+1:o.Nnet)'];

% initialize mu,sdev. for PCA net, we do not output any of them.
muChan1 = []; sdevChan1 = []; scale = []; offset = []; outputSize = [];
netFilesFolder = './networks/NEWn256Mtimes24modesFFTNets/velPredPCAin_mode';

for imode = 1 : o.nVelModes
  pmode = o.velActiveModes(imode);
  load([netFilesFolder num2str(pmode) '_net_w1step.mat'])
  FCnets{imode} = net;
end

o.MVnets = FCnets; o.muChan_MV = muChan1; o.sdevChan_MV = sdevChan1;
o.scale_MV = scale; o.offset_MV = offset; o.MVoutSize = outputSize;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = loadTenNetFiles(o)
% network for inverse of tension matrix on self bending
o.tenPredModes = [(1:o.nTenModes/2)';(o.Nnet-o.nTenModes/2+1:o.Nnet)'];
load('./networks/NEWfcPCAtenMaTonBendN256_32modes_FCNet_w1step.mat')
o.tenBendNets = net; o.muChan_tenBend = muChan1; 
o.sdevChan_tenBend = sdevChan1; o.scale_tenBend = scale; 
o.offset_tenBend = offset; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = loadTenVnetFiles(o)
% network for inverse of tension matrix on vback

% Fourier modes
o.TenVactiveModes = [(1:o.nVelModes/2)';(o.Nnet-o.nVelModes/2+1:o.Nnet)'];

% initialize mu,sdev. for PCA net, we do not output any of them.
muChan1 = []; sdevChan1 = []; scale = []; offset = []; outputSize = [];
netFilesFolder = './networks/NEWn256tenMatTimes24modesFFTNets/velPredPCAin_mode';

for imode = 1 : o.nVelModes
  pmode = o.TenVactiveModes(imode);
  load([netFilesFolder num2str(pmode) '_net_w1step.mat'])
  FCnets{imode} = net;
end

o.tenVnets = FCnets; o.muChan_tenV = muChan1; o.sdevChan_tenV = sdevChan1;
o.scale_tenV = scale; o.offset_tenV = offset; o.tenVoutSize = outputSize;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = initWalls(o,Xwalls,XwallsInt,XwallsExt)
tt = o.tt;

if ~isempty(Xwalls) % if not DLD
nvbd = numel(Xwalls(1,:)); % number of walls
opWall = tt.opWall;
% velocity on solid walls
[uwalls,~] = tt.farField(Xwalls,[]);
% build walls
walls = capsules(Xwalls,[],uwalls,zeros(nvbd,1),zeros(nvbd,1),1);
walls.setUpRate();
o.walls = walls;
if nvbd > 1
  o.NearW2W = walls.getZone([],1);
end

% build the double-layer potential matrix for walls and save on memory
tt.wallDLP = opWall.stokesDLmatrix(walls);
tt.wallDLPnoCorr = opWall.stokesDLmatrixNoCorr(walls);

% N0 to remove rank-1 deficiency
tt.wallN0 = opWall.stokesN0matrix(walls);

% inverse of the walls matrix (-1/2*density + DLP + N0 + R-S)
tt.matFreeWalls = false; % so, create walls matrix and save it
tt.bdiagWall = tt.wallsPrecond(walls); % inverse of the walls matrix

% assign empty for DLD related objects
o.wallsInt = []; o.wallsExt = [];
o.NearWint2Wint = [];
o.NearWext2Wint = [];
o.NearWint2Wext = [];
else % then DLD 
% poten classes for walls
potWallInt = tt.opWallInt;
potWallExt = tt.opWallExt;
% velocity on solid walls coming from no slip boundary condition
[uwallsExt,uwallsInt] = tt.farField(XwallsExt,XwallsInt);

% number of pillars
nvbdInt = numel(XwallsInt(1,:));
% build walls
wallsInt = capsules(XwallsInt,[],uwallsInt,...
      zeros(nvbdInt,1),zeros(nvbdInt,1),1);
wallsExt = capsules(XwallsExt,[],uwallsExt,0,0,1);
wallsInt.setUpRate(potWallInt);
wallsExt.setUpRate(potWallExt);
o.wallsInt = wallsInt; o.wallsExt = wallsExt;

% get the near structure
o.NearWint2Wint = wallsInt.getZone([],1);
o.NearWext2Wint = wallsExt.getZone(wallsInt,2);
o.NearWint2Wext = wallsInt.getZone(wallsExt,2);

% build the double-layer potential matrix for walls and save on memory
% build the double layer potential matrix for walls and save on memory
if isempty(tt.wallDLPint)
  tt.wallDLPint = potWallInt.stokesDLmatrix(wallsInt);
end
if isempty(tt.wallDLPext)
  tt.wallDLPext = potWallExt.stokesDLmatrix(wallsExt);
end

if isempty(tt.wallDLPintNoCorr)
  tt.wallDLPintNoCorr = potWallInt.stokesDLmatrixNoCorr(wallsInt);
end
if isempty(tt.wallDLPextNoCorr)
  tt.wallDLPextNoCorr = potWallExt.stokesDLmatrixNoCorr(wallsExt);    
end


% N0 to remove rank-1 deficiency
if isempty(tt.wallN0)
  tt.wallN0 = potWallExt.stokesN0matrix(wallsExt);
end

% inverse of the walls matrix (-1/2*density + DLP + N0 + R-S)
tt.matFreeWalls = false; % so, create walls matrix and save it
tt.bdiagWall = tt.wallsPrecondDiffDisc(wallsInt,wallsExt);

% assign empty for walls related objects
o.walls = [];
o.NearW2W = [];
end % if 


o.tt = tt; 
% tt includes tt.wallDLPandRSmat which is walls matrix saved and inverse of walls matrix

end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xnew,tenNew,etaNew,etaIntNew,etaExtNew,RSnew] = DNNsolve(o,...
        Xold,tenOld,etaOld,etaIntOld,etaExtOld,RSold)

tt = o.tt;
op = tt.op;
if o.confined
  if isempty(o.wallsInt) % not DLD
  opWall = tt.opWall;
  walls = o.walls;
  uWalls = walls.u;
  nvbd = walls.nv;
  Nbd = walls.N;
  diffDiscWalls = false;
  etaIntNew = []; etaExtNew = [];
  else
  opWallInt = tt.opWallInt;
  opWallExt = tt.opWallExt;  
  wallsInt = o.wallsInt;
  wallsExt = o.wallsExt;
  nvbdInt = wallsInt.nv;
  nvbdExt = 1;
  nvbd = nvbdExt + nvbdInt;
  NbdInt = wallsInt.N;
  NbdExt = wallsExt.N;
  diffDiscWalls = true;
  etaNew = [];
  end % if isempty(o.wallsInt)
  RSnew = zeros(3,nvbd);
else
  vback = o.vinf(Xold);
  etaNew = []; RSnew = []; etaIntNew = []; etaExtNew = [];
end

% explicit time stepping w/ splitting
vesicle = capsules(Xold,[],[],o.kappa,1,1);
vesicle.setUpRate();
nv = vesicle.nv;
N = vesicle.N;

if tt.fmm
  Nup = op.LPuprate * vesicle.N;
  Xup = [interpft(Xold(1:end/2,:),Nup);interpft(Xold(end/2+1:end,:),Nup)];
  vesicleUp = capsules(Xup,[],[],vesicle.kappa,vesicle.viscCont,0); 
  SLPnoCorr = op.stokesSLmatrixNoCorr(vesicleUp);
else
  SLPnoCorr = [];
end

% Get the near structure
if o.confined
  if ~diffDiscWalls  
    [NearV2V,NearV2W] = vesicle.getZone(walls,3);
    [~,NearW2V] = walls.getZone(vesicle,2);
    if nvbd > 1
      NearW2W = o.NearW2W;  
    end % if nvbd > 1
  else
    % Get the near structure, so that we know what to neglect
    [NearV2V,NearV2Wint] = vesicle.getZone(wallsInt,3);
      
    if isempty(tt.NearWint2Wint)
      [tt.NearWint2Wint,NearWint2V] = wallsInt.getZone(vesicle,3);
    else
      % there is no need to compute W2W again, since they do not move
      [~,NearWint2V] = wallsInt.getZone(vesicle,2);
    end
    NearWint2Wint = tt.NearWint2Wint;

    [~,NearV2Wext] = vesicle.getZone(wallsExt,2);
    [~,NearWext2V] = wallsExt.getZone(vesicle,2);

    if isempty(tt.NearWint2Wext)
      % there is no need to compute W2W again, since they do not move
      [~,tt.NearWint2Wext] = wallsInt.getZone(wallsExt,2);
    end
    NearWint2Wext = tt.NearWint2Wext;

    if isempty(tt.NearWext2Wint)
      % there is no need to compute W2W again, since they do not move
      [~,tt.NearWext2Wint] = wallsExt.getZone(wallsInt,2);    
    end
    NearWext2Wint = tt.NearWext2Wint;
  end % if ~diffDiscWalls
else
  NearV2V = vesicle.getZone([],1);    
end % if o.confined


% 1) EXPLICIT TENSION AT THE CURRENT STEP

% Compute bending forces + old tension forces
fBend = vesicle.tracJump(Xold,zeros(N,nv));
fTen = vesicle.tracJump(zeros(2*N,nv),tenOld);
tracJump = fBend+fTen;

if nv > 1
  % SLP due to traction jump
  kernelRegul = @op.exactStokesSLregul; % regularized kernel for near-field
  kernelDirect = @op.exactStokesSL; % direct evaluation for far-field
  if tt.fmm
    kernel = @op.exactStokesSLfmm;
  else
    kernel = @op.exactStokesSLregul;
  end
  ves2vesSLP = op.evaluateLayerPoten(vesicle,tracJump,SLPnoCorr,NearV2V,...
        kernel,kernelDirect,kernelRegul,vesicle,true);
else
  ves2vesSLP = zeros(2*N,1);
end

% DLP due to walls or background flow
if o.confined
  if ~diffDiscWalls  
    kernelRegul = @opWall.exactStokesDLregul;
    kernelDirect = @opWall.exactStokesDL;
    if tt.fmmDLP
      kernel = @opWall.exactStokesDLnewfmm;
    else
      kernel = @opWall.exactStokesDLregul;
    end
    vback = opWall.evaluateLayerPoten(walls,etaOld,[],NearW2V,...
        kernel,kernelDirect,kernelRegul,vesicle,false);  
  else % then DLD
    kernelRegul = @opWallInt.exactStokesDLregul;
    kernelRegul2 = @opWallExt.exactStokesDLregul;
    kernelDirect = @opWallInt.exactStokesDL;
    kernelDirect2 = @opWallExt.exactStokesDL;
    if tt.fmmDLP
      kernel = @opWallInt.exactStokesDLnewfmm;
      kernel2 = @opWallExt.exactStokesDLnewfmm;
    else
      kernel = @opWallInt.exactStokesDLregul;
      kernel2 = @opWallExt.exactStokesDLregul;
    end
    % interior walls to vesicles
    wallInt2vesDLP = opWallInt.evaluateLayerPoten(wallsInt,etaIntOld,[],NearWint2V,...
        kernel,kernelDirect,kernelRegul,vesicle,false);
    % exterior wall to vesicles
    wallExt2vesDLP = opWallExt.evaluateLayerPoten(wallsExt,etaExtOld,[],NearWext2V,...
        kernel2,kernelDirect2,kernelRegul2,vesicle,false);
    
    vback = wallInt2vesDLP + wallExt2vesDLP;
  end %if ~diffDiscWalls
  
  for k = 2:nvbd
    stokeslet = RSold(1:2,k); rotlet = RSold(3,k);  
    if ~diffDiscWalls
      vback = vback + tt.RSlets(vesicle.X,walls.center(:,k),stokeslet,rotlet);
    else
      vback = vback + tt.RSlets(vesicle.X,wallsInt.center(:,k-1),stokeslet,rotlet);
    end
  end % k = 2:nvbd (stokeslet and rotlet)
  
end % if o.confined

scaling = zeros(nv,1); rotate = zeros(nv,1); trans = zeros(2,nv);
sortIdx = zeros(o.Nnet,nv);
XinputFourier = zeros(o.nComp,1,1,nv);
XinputTenBen = zeros(o.nComp,1,1,nv);
for k = 1 : nv
  % Standardize vesicle  
  [Xstand,scaling(k),rotate(k),trans(:,k),sortIdx(:,k)] = ...
      o.standardizationStep(Xold(:,k));
  % Prepare inputs for Fourier network and self-bending network
  XinputFourier(:,:,:,k) = o.prepareInputForNet(Xstand,'tensionOnFourier');
  XinputTenBen(:,:,:,k) = o.prepareInputForNet(Xstand,'tensionSelfBend');
end
% Approximate inv(Div*G*Ten)*Div*vExt
vBackSolve = o.invTenMatOnVback(XinputFourier,vback+ves2vesSLP,...
    scaling,rotate,sortIdx);
% Approximate inv(Div*G*Ten)*G*(-Ben)*x
selfBendSolve = o.invTenMatOnSelfBend(XinputTenBen,scaling,sortIdx,N);
% update tension
tenNew = -(vBackSolve+selfBendSolve);

% Update the traction jump
fTen = vesicle.tracJump(zeros(2*N,nv),tenNew); tracJump = fBend+fTen;

% 2) SOLVE FOR DENSITY and RS ON WALLS and UPDATE wall2vesDLP
if o.confined
  % vesicle2wall interactions
  kernelRegul = @op.exactStokesSLregul;
  kernelDirect = @op.exactStokesSL;
  if tt.fmm
    kernel = @op.exactStokesSLfmm;
  else
    kernel = @op.exactStokesSLregul;
  end
  
  % SLP due to vesicles on walls
  if ~diffDiscWalls
    ves2wallSLP = op.evaluateLayerPoten(vesicle,tracJump,[],NearV2W,...
        kernel,kernelDirect,kernelRegul,walls,0);
  
    RHS = [uWalls(:)-ves2wallSLP(:);zeros(3*(nvbd-1),1)];
    etaRS = tt.bdiagWall*RHS;
    etaNew = zeros(2*Nbd,nvbd);
    for iw = 1 : nvbd
      etaNew(:,iw) = etaRS((iw-1)*2*Nbd+1:iw*2*Nbd);  
      if iw <= nvbd-1
        RSnew(:,iw+1) = etaRS(2*Nbd*nvbd+(iw-1)*3+1:2*Nbd*nvbd+iw*3);
      end
    end 
  else % if DLD  
    ves2wallInt = op.evaluateLayerPoten(vesicle,tracJump,[],NearV2Wint,...
        kernel,kernelDirect,kernelRegul,wallsInt,0);
    ves2wallExt = op.evaluateLayerPoten(vesicle,tracJump,[],NearV2Wext,...
        kernel,kernelDirect,kernelRegul,wallsExt,0);
    
    RHSint = wallsInt.u(:)-ves2wallInt(:);
    RHSext = wallsExt.u(:)-ves2wallExt(:);
    RHS = [RHSext;RHSint;zeros(3*(nvbd-1),1)];
    etaRS = tt.bdiagWall*RHS;
    etaExtNew = etaRS(1:2*NbdExt);
    etaRS = etaRS(2*NbdExt+1:end);
    etaIntNew = zeros(2*NbdInt,nvbdInt);
    for iw = 1 : nvbdInt
      etaIntNew(:,iw) = etaRS((iw-1)*2*NbdInt+1:iw*2*NbdInt);  
      if iw <= nvbdInt
        RSnew(:,iw+1) = etaRS(2*NbdInt*nvbdInt+(iw-1)*3+1:2*NbdInt*nvbdInt+iw*3);
      end
    end   
  end % if ~diffDiscWalls
  
  % Update vback due to new density,rotlets and stokeslet
  if ~diffDiscWalls  
    kernelRegul = @opWall.exactStokesDLregul;
    kernelDirect = @opWall.exactStokesDL;
    if tt.fmmDLP
      kernel = @opWall.exactStokesDLnewfmm;
    else
      kernel = @opWall.exactStokesDLregul;
    end
    vback = opWall.evaluateLayerPoten(walls,etaNew,[],NearW2V,...
        kernel,kernelDirect,kernelRegul,vesicle,false);  
  else % then DLD
    kernelRegul = @opWallInt.exactStokesDLregul;
    kernelRegul2 = @opWallExt.exactStokesDLregul;
    kernelDirect = @opWallInt.exactStokesDL;
    kernelDirect2 = @opWallExt.exactStokesDL;
    if tt.fmmDLP
      kernel = @opWallInt.exactStokesDLnewfmm;
      kernel2 = @opWallExt.exactStokesDLnewfmm;
    else
      kernel = @opWallInt.exactStokesDLregul;
      kernel2 = @opWallExt.exactStokesDLregul;
    end
    % interior walls to vesicles
    wallInt2vesDLP = opWallInt.evaluateLayerPoten(wallsInt,etaIntNew,[],NearWint2V,...
        kernel,kernelDirect,kernelRegul,vesicle,false);
    % exterior wall to vesicles
    wallExt2vesDLP = opWallExt.evaluateLayerPoten(wallsExt,etaExtNew,[],NearWext2V,...
        kernel2,kernelDirect2,kernelRegul2,vesicle,false);
    
    vback = wallInt2vesDLP + wallExt2vesDLP;
  end %if ~diffDiscWalls
  
  for k = 2:nvbd
    stokeslet = RSnew(1:2,k); rotlet = RSnew(3,k);  
    if ~diffDiscWalls
      vback = vback + tt.RSlets(vesicle.X,walls.center(:,k),stokeslet,rotlet);
    else
      vback = vback + tt.RSlets(vesicle.X,wallsInt.center(:,k-1),stokeslet,rotlet);
    end
  end % k = 2:nvbd (stokeslet and rotlet)
end % if iconfined    

% 3) SOLVE FOR NEXT TIME STEP USING OPERATOR SPLITTING
% EXPLICIT-TREATMENT OF INTER-VESICLE BENDING AND TENSION
% TENSION AND X are coupled (semi-implicit).
if nv > 1
  % SLP due to traction jump
  kernelRegul = @op.exactStokesSLregul; % regularized kernel for near-field
  kernelDirect = @op.exactStokesSL; % direct evaluation for far-field
  if tt.fmm
    kernel = @op.exactStokesSLfmm;
  else
    kernel = @op.exactStokesSLregul;
  end
  ves2vesSLP = op.evaluateLayerPoten(vesicle,tracJump,SLPnoCorr,NearV2V,...
        kernel,kernelDirect,kernelRegul,vesicle,true);
else
  ves2vesSLP = zeros(2*N,1);
end

% Walls2Ves is already computed if confined (vback)    
rotate = zeros(nv,1); trans = zeros(2,nv); sortIdx = zeros(o.Nnet,nv);
scaling = zeros(nv,1); Xinput = zeros(o.nComp,1,1,nv);
for k = 1 : nv
  % 1) TRANSLATION W/ NETWORK
  % Standardize vesicle
  [Xstand,scaling(k),rotate(k),trans(:,k),sortIdx(:,k)] = ...
      o.standardizationStep(Xold(:,k));
  % Prepare input for advection network
  Xinput(:,1,1,k) = o.prepareInputForNet(Xstand,'advection');
end
% Take a step due to advaction
Xmid = o.translateVinfwNN(Xinput,vback + ves2vesSLP,Xold,rotate,sortIdx);

if ~o.variableKbDt
  rotate = zeros(nv,1); trans = zeros(2,nv); sortIdx = zeros(o.Nnet,nv);
  scaling = zeros(nv,1); Xinput = zeros(o.nCompRelax,1,1,nv);
  for k = 1 : nv
    % 2) RELAXATION w/ NETWORK
    % Standardize vesicle Xmid
    [Xstand,scaling(k),rotate(k),trans(:,k),sortIdx(:,k)] = ...
        o.standardizationStep(Xmid(:,k));
    % Prepare input for relaxation network
    Xinput(:,1,1,k) = o.prepareInputForNet(Xstand,'relaxation');
  end
  % Take a step
  Xnew = o.relaxWNN(Xinput,scaling,rotate,trans,sortIdx,N);
else
  Xnew = o.relaxWNNvariableKbDt(Xmid);    
end

end % DNNsolve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xnew,tenNew,etaNew,etaIntNew,etaExtNew,RSnew] = DNNsolveAltern(o,...
        Xold,tenOld,etaOld,etaIntOld,etaExtOld,RSold)

tt = o.tt;
op = tt.op;
if o.confined
  if isempty(o.wallsInt) % not DLD
  opWall = tt.opWall;
  walls = o.walls;
  uWalls = walls.u;
  nvbd = walls.nv;
  Nbd = walls.N;
  diffDiscWalls = false;
  etaIntNew = []; etaExtNew = [];
  else
  opWallInt = tt.opWallInt;
  opWallExt = tt.opWallExt;  
  wallsInt = o.wallsInt;
  wallsExt = o.wallsExt;
  nvbdInt = wallsInt.nv;
  nvbdExt = 1;
  nvbd = nvbdExt + nvbdInt;
  NbdInt = wallsInt.N;
  NbdExt = wallsExt.N;
  diffDiscWalls = true;
  etaNew = [];
  end % if isempty(o.wallsInt)
  RSnew = zeros(3,nvbd);
else
  vback = o.vinf(Xold);
  etaNew = []; RSnew = []; etaIntNew = []; etaExtNew = [];
end

% explicit time stepping w/ splitting
vesicle = capsules(Xold,[],[],o.kappa,1,1);
vesicle.setUpRate();
nv = vesicle.nv;
N = vesicle.N;

if tt.fmm
  Nup = op.LPuprate * vesicle.N;
  Xup = [interpft(Xold(1:end/2,:),Nup);interpft(Xold(end/2+1:end,:),Nup)];
  vesicleUp = capsules(Xup,[],[],vesicle.kappa,vesicle.viscCont,0); 
  SLPnoCorr = op.stokesSLmatrixNoCorr(vesicleUp);
else
  SLPnoCorr = [];
end

% Get the near structure
if o.confined
  if ~diffDiscWalls  
    [NearV2V,NearV2W] = vesicle.getZone(walls,3);
    [~,NearW2V] = walls.getZone(vesicle,2);
    if nvbd > 1
      NearW2W = o.NearW2W;  
    end % if nvbd > 1
  else
    % Get the near structure, so that we know what to neglect
    [NearV2V,NearV2Wint] = vesicle.getZone(wallsInt,3);
      
    if isempty(tt.NearWint2Wint)
      [tt.NearWint2Wint,NearWint2V] = wallsInt.getZone(vesicle,3);
    else
      % there is no need to compute W2W again, since they do not move
      [~,NearWint2V] = wallsInt.getZone(vesicle,2);
    end
    NearWint2Wint = tt.NearWint2Wint;

    [~,NearV2Wext] = vesicle.getZone(wallsExt,2);
    [~,NearWext2V] = wallsExt.getZone(vesicle,2);

    if isempty(tt.NearWint2Wext)
      % there is no need to compute W2W again, since they do not move
      [~,tt.NearWint2Wext] = wallsInt.getZone(wallsExt,2);
    end
    NearWint2Wext = tt.NearWint2Wext;

    if isempty(tt.NearWext2Wint)
      % there is no need to compute W2W again, since they do not move
      [~,tt.NearWext2Wint] = wallsExt.getZone(wallsInt,2);    
    end
    NearWext2Wint = tt.NearWext2Wint;
  end % if ~diffDiscWalls
else
  NearV2V = vesicle.getZone([],1);    
end % if o.confined


% 1) SOLVE FOR VESICLES' NEW POSITION WITH OLD TENSION AND DENSITY

% Compute bending forces + old tension forces
fBend = vesicle.tracJump(Xold,zeros(N,nv));
fTen = vesicle.tracJump(zeros(2*N,nv),tenOld);
tracJump = fBend+fTen;

if nv > 1
  % SLP due to traction jump
  kernelRegul = @op.exactStokesSLregul; % regularized kernel for near-field
  kernelDirect = @op.exactStokesSL; % direct evaluation for far-field
  if tt.fmm
    kernel = @op.exactStokesSLfmm;
  else
    kernel = @op.exactStokesSLregul;
  end
  ves2vesSLP = op.evaluateLayerPoten(vesicle,tracJump,SLPnoCorr,NearV2V,...
        kernel,kernelDirect,kernelRegul,vesicle,true);
else
  ves2vesSLP = zeros(2*N,1);
end

% DLP due to walls or background flow
if o.confined
  if ~diffDiscWalls  
    kernelRegul = @opWall.exactStokesDLregul;
    kernelDirect = @opWall.exactStokesDL;
    if tt.fmmDLP
      kernel = @opWall.exactStokesDLnewfmm;
    else
      kernel = @opWall.exactStokesDLregul;
    end
    vback = opWall.evaluateLayerPoten(walls,etaOld,[],NearW2V,...
        kernel,kernelDirect,kernelRegul,vesicle,false);  
  else % then DLD
    kernelRegul = @opWallInt.exactStokesDLregul;
    kernelRegul2 = @opWallExt.exactStokesDLregul;
    kernelDirect = @opWallInt.exactStokesDL;
    kernelDirect2 = @opWallExt.exactStokesDL;
    if tt.fmmDLP
      kernel = @opWallInt.exactStokesDLnewfmm;
      kernel2 = @opWallExt.exactStokesDLnewfmm;
    else
      kernel = @opWallInt.exactStokesDLregul;
      kernel2 = @opWallExt.exactStokesDLregul;
    end
    % interior walls to vesicles
    wallInt2vesDLP = opWallInt.evaluateLayerPoten(wallsInt,etaIntOld,[],NearWint2V,...
        kernel,kernelDirect,kernelRegul,vesicle,false);
    % exterior wall to vesicles
    wallExt2vesDLP = opWallExt.evaluateLayerPoten(wallsExt,etaExtOld,[],NearWext2V,...
        kernel2,kernelDirect2,kernelRegul2,vesicle,false);
    
    vback = wallInt2vesDLP + wallExt2vesDLP;
  end %if ~diffDiscWalls
  
  for k = 2:nvbd
    stokeslet = RSold(1:2,k); rotlet = RSold(3,k);  
    if ~diffDiscWalls
      vback = vback + tt.RSlets(vesicle.X,walls.center(:,k),stokeslet,rotlet);
    else
      vback = vback + tt.RSlets(vesicle.X,wallsInt.center(:,k-1),stokeslet,rotlet);
    end
  end % k = 2:nvbd (stokeslet and rotlet)
  
end % if o.confined

rotate = zeros(nv,1); trans = zeros(2,nv); sortIdx = zeros(o.Nnet,nv);
scaling = zeros(nv,1); Xinput = zeros(o.nComp,1,1,nv);
for k = 1 : nv
  % 1) TRANSLATION W/ NETWORK
  % Standardize vesicle
  [Xstand,scaling(k),rotate(k),trans(:,k),sortIdx(:,k)] = ...
      o.standardizationStep(Xold(:,k));
  % Prepare input for advection network
  Xinput(:,1,1,k) = o.prepareInputForNet(Xstand,'advection');
end
% Take a step due to advaction
Xmid = o.translateVinfwNN(Xinput,vback + ves2vesSLP,Xold,rotate,sortIdx);

if ~o.variableKbDt
  rotate = zeros(nv,1); trans = zeros(2,nv); sortIdx = zeros(o.Nnet,nv);
  scaling = zeros(nv,1); Xinput = zeros(o.nCompRelax,1,1,nv);
  for k = 1 : nv
    % 2) RELAXATION w/ NETWORK
    % Standardize vesicle Xmid
    [Xstand,scaling(k),rotate(k),trans(:,k),sortIdx(:,k)] = ...
        o.standardizationStep(Xmid(:,k));
    % Prepare input for relaxation network
    Xinput(:,1,1,k) = o.prepareInputForNet(Xstand,'relaxation');
  end
  % Take a step
  Xnew = o.relaxWNN(Xinput,scaling,rotate,trans,sortIdx,N);
else
  Xnew = o.relaxWNNvariableKbDt(Xmid);    
end

% 2) EXPLICIT TENSION AT THE CURRENT STEP
vesicle = capsules(Xnew,[],[],o.kappa,1,1);
vesicle.setUpRate();

if tt.fmm
  Nup = op.LPuprate * vesicle.N;
  Xup = [interpft(Xnew(1:end/2,:),Nup);interpft(Xnew(end/2+1:end,:),Nup)];
  vesicleUp = capsules(Xup,[],[],vesicle.kappa,vesicle.viscCont,0); 
  SLPnoCorr = op.stokesSLmatrixNoCorr(vesicleUp);
else
  SLPnoCorr = [];
end

% Get the near structure
if o.confined
  if ~diffDiscWalls  
    [NearV2V,NearV2W] = vesicle.getZone(walls,3);
    [~,NearW2V] = walls.getZone(vesicle,2);
    if nvbd > 1
      NearW2W = o.NearW2W;  
    end % if nvbd > 1
  else
    % Get the near structure, so that we know what to neglect
    [NearV2V,NearV2Wint] = vesicle.getZone(wallsInt,3);
      
    if isempty(tt.NearWint2Wint)
      [tt.NearWint2Wint,NearWint2V] = wallsInt.getZone(vesicle,3);
    else
      % there is no need to compute W2W again, since they do not move
      [~,NearWint2V] = wallsInt.getZone(vesicle,2);
    end
    NearWint2Wint = tt.NearWint2Wint;

    [~,NearV2Wext] = vesicle.getZone(wallsExt,2);
    [~,NearWext2V] = wallsExt.getZone(vesicle,2);

    if isempty(tt.NearWint2Wext)
      % there is no need to compute W2W again, since they do not move
      [~,tt.NearWint2Wext] = wallsInt.getZone(wallsExt,2);
    end
    NearWint2Wext = tt.NearWint2Wext;

    if isempty(tt.NearWext2Wint)
      % there is no need to compute W2W again, since they do not move
      [~,tt.NearWext2Wint] = wallsExt.getZone(wallsInt,2);    
    end
    NearWext2Wint = tt.NearWext2Wint;
  end % if ~diffDiscWalls
else
  NearV2V = vesicle.getZone([],1);    
end % if o.confined

% Compute bending forces + old tension forces
fBend = vesicle.tracJump(Xnew,zeros(N,nv));
tracJump = fBend+fTen;

if nv > 1
  % SLP due to traction jump
  kernelRegul = @op.exactStokesSLregul; % regularized kernel for near-field
  kernelDirect = @op.exactStokesSL; % direct evaluation for far-field
  if tt.fmm
    kernel = @op.exactStokesSLfmm;
  else
    kernel = @op.exactStokesSLregul;
  end
  ves2vesSLP = op.evaluateLayerPoten(vesicle,tracJump,SLPnoCorr,NearV2V,...
        kernel,kernelDirect,kernelRegul,vesicle,true);
else
  ves2vesSLP = zeros(2*N,1);
end

% DLP due to walls or background flow
if o.confined
  if ~diffDiscWalls  
    kernelRegul = @opWall.exactStokesDLregul;
    kernelDirect = @opWall.exactStokesDL;
    if tt.fmmDLP
      kernel = @opWall.exactStokesDLnewfmm;
    else
      kernel = @opWall.exactStokesDLregul;
    end
    vback = opWall.evaluateLayerPoten(walls,etaOld,[],NearW2V,...
        kernel,kernelDirect,kernelRegul,vesicle,false);  
  else % then DLD
    kernelRegul = @opWallInt.exactStokesDLregul;
    kernelRegul2 = @opWallExt.exactStokesDLregul;
    kernelDirect = @opWallInt.exactStokesDL;
    kernelDirect2 = @opWallExt.exactStokesDL;
    if tt.fmmDLP
      kernel = @opWallInt.exactStokesDLnewfmm;
      kernel2 = @opWallExt.exactStokesDLnewfmm;
    else
      kernel = @opWallInt.exactStokesDLregul;
      kernel2 = @opWallExt.exactStokesDLregul;
    end
    % interior walls to vesicles
    wallInt2vesDLP = opWallInt.evaluateLayerPoten(wallsInt,etaIntOld,[],NearWint2V,...
        kernel,kernelDirect,kernelRegul,vesicle,false);
    % exterior wall to vesicles
    wallExt2vesDLP = opWallExt.evaluateLayerPoten(wallsExt,etaExtOld,[],NearWext2V,...
        kernel2,kernelDirect2,kernelRegul2,vesicle,false);
    
    vback = wallInt2vesDLP + wallExt2vesDLP;
  end %if ~diffDiscWalls
  
  for k = 2:nvbd
    stokeslet = RSold(1:2,k); rotlet = RSold(3,k);  
    if ~diffDiscWalls
      vback = vback + tt.RSlets(vesicle.X,walls.center(:,k),stokeslet,rotlet);
    else
      vback = vback + tt.RSlets(vesicle.X,wallsInt.center(:,k-1),stokeslet,rotlet);
    end
  end % k = 2:nvbd (stokeslet and rotlet)
  
end % if o.confined

scaling = zeros(nv,1); rotate = zeros(nv,1); trans = zeros(2,nv);
sortIdx = zeros(o.Nnet,nv);
XinputFourier = zeros(o.nComp,1,1,nv);
XinputTenBen = zeros(o.nComp,1,1,nv);
for k = 1 : nv
  % Standardize vesicle  
  [Xstand,scaling(k),rotate(k),trans(:,k),sortIdx(:,k)] = ...
      o.standardizationStep(Xnew(:,k));
  % Prepare inputs for Fourier network and self-bending network
  XinputFourier(:,:,:,k) = o.prepareInputForNet(Xstand,'tensionOnFourier');
  XinputTenBen(:,:,:,k) = o.prepareInputForNet(Xstand,'tensionSelfBend');
end
% Approximate inv(Div*G*Ten)*Div*vExt
vBackSolve = o.invTenMatOnVback(XinputFourier,vback+ves2vesSLP,...
    scaling,rotate,sortIdx);
% Approximate inv(Div*G*Ten)*G*(-Ben)*x
selfBendSolve = o.invTenMatOnSelfBend(XinputTenBen,scaling,sortIdx,N);
% update tension
tenNew = -(vBackSolve+selfBendSolve);

% 3) SOLVE FOR DENSITY and RS ON WALLS and UPDATE wall2vesDLP
if o.confined
  % Update the traction jump
  fTen = vesicle.tracJump(zeros(2*N,nv),tenNew); tracJump = fBend+fTen;  
  
  % vesicle2wall interactions
  kernelRegul = @op.exactStokesSLregul;
  kernelDirect = @op.exactStokesSL;
  if tt.fmm
    kernel = @op.exactStokesSLfmm;
  else
    kernel = @op.exactStokesSLregul;
  end
  
  % SLP due to vesicles on walls
  if ~diffDiscWalls
    ves2wallSLP = op.evaluateLayerPoten(vesicle,tracJump,[],NearV2W,...
        kernel,kernelDirect,kernelRegul,walls,0);
  
    RHS = [uWalls(:)-ves2wallSLP(:);zeros(3*(nvbd-1),1)];
    etaRS = tt.bdiagWall*RHS;
    etaNew = zeros(2*Nbd,nvbd);
    for iw = 1 : nvbd
      etaNew(:,iw) = etaRS((iw-1)*2*Nbd+1:iw*2*Nbd);  
      if iw <= nvbd-1
        RSnew(:,iw+1) = etaRS(2*Nbd*nvbd+(iw-1)*3+1:2*Nbd*nvbd+iw*3);
      end
    end 
  else % if DLD  
    ves2wallInt = op.evaluateLayerPoten(vesicle,tracJump,[],NearV2Wint,...
        kernel,kernelDirect,kernelRegul,wallsInt,0);
    ves2wallExt = op.evaluateLayerPoten(vesicle,tracJump,[],NearV2Wext,...
        kernel,kernelDirect,kernelRegul,wallsExt,0);
    
    RHSint = wallsInt.u(:)-ves2wallInt(:);
    RHSext = wallsExt.u(:)-ves2wallExt(:);
    RHS = [RHSext;RHSint;zeros(3*(nvbd-1),1)];
    etaRS = tt.bdiagWall*RHS;
    etaExtNew = etaRS(1:2*NbdExt);
    etaRS = etaRS(2*NbdExt+1:end);
    etaIntNew = zeros(2*NbdInt,nvbdInt);
    for iw = 1 : nvbdInt
      etaIntNew(:,iw) = etaRS((iw-1)*2*NbdInt+1:iw*2*NbdInt);  
      if iw <= nvbdInt
        RSnew(:,iw+1) = etaRS(2*NbdInt*nvbdInt+(iw-1)*3+1:2*NbdInt*nvbdInt+iw*3);
      end
    end   
  end % if ~diffDiscWalls
end % if iconfined    

end % DNNsolveAltern
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vBackSolve = invTenMatOnVback(o,Xinput,vext,scaling,rotate,sortIdx)
% Approximate inv(Div*G*Ten)*Div*vExt 
    
% number of vesicles
nv = numel(vext(1,:));
% number of points of exact solve
N = numel(vext(:,1))/2;    
% number of points for network
Nnet = o.Nnet;
activeModes = o.TenVactiveModes;
outputSize = o.tenVoutSize;
nets = o.tenVnets;

% Upsample vext
vext = [interpft(vext(1:end/2,:),Nnet);interpft(vext(end/2+1:end,:),Nnet)];

% Initialize the multiplication matrix Z = inv(DivGT)DivPhi_k
Z1 = zeros(Nnet,numel(activeModes),nv); Z2 = Z1;

for k = 1 : numel(activeModes)
  pred = predict(nets{k},Xinput)';
  Z1(:,k,:) = interpft(pred(1:outputSize/2,:),Nnet);
  Z2(:,k,:) = interpft(pred(outputSize/2+1:outputSize,:),Nnet);
end

vBackSolve = zeros(N,nv);
for k = 1 : nv
  % Take fft of the velocity, standardize velocity
  vextStand = o.standardize(vext(:,k),[0;0],rotate(k),1,sortIdx(:,k));
  z = vextStand(1:end/2)+1i*vextStand(end/2+1:end);
  zh = fft(z);
  V1 = real(zh(activeModes)); V2 = imag(zh(activeModes));

  % Compute the approximation to inv(Div*G*Ten)*Div*vExt
  MVextStand = (Z1(:,:,k)*V1+Z2(:,:,k)*V2);

  % Destandardize the multiplication
  MVext = zeros(size(MVextStand));
  MVext(sortIdx(:,k)) = MVextStand;

  % downsample MVext
  vBackSolve(:,k) = interpft(MVext,N);
end
end % invTenMatOnVback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function selfBendSolve = invTenMatOnSelfBend(o,Xinput,scaling,sortIdx,N)
% Approximate inv(Div*G*Ten)*G*(-Ben)*x
% number of vesicles
nv = numel(sortIdx(1,:));
% number of points for network
Nnet = o.Nnet;
% network and active Fourier modes to be predicted
net = o.tenBendNets;
activeModes = o.tenPredModes;

% predict the active modes
realImagActive = predict(net,Xinput)';
% initialize full output
zFull = zeros(Nnet,nv);
zFull(activeModes,:) = realImagActive(1:end/2,:)+1i*realImagActive(end/2+1:end,:);

% ifft and find upsampled output
upsampOutStand = real(ifft(zFull)*Nnet);

% destandardize
selfBendSolve = zeros(N,nv);
for k = 1 : nv
  upsampOut = zeros(Nnet,1);  
  upsampOut(sortIdx(:,k)) = o.kappa*upsampOutStand(:,k)/scaling(k)^2;
  % downsample
  selfBendSolve(:,k) = interpft(upsampOut,N);
end

end % invTenMatOnSelfBend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xnew = relaxWNN(o,Xinput,scaling,rotate,trans,sortIdx,N)
% load network
nets = o.bendNets; muChan1 = o.muChan_bend; sdevChan1 = o.sdevChan_bend;
scale = o.scale_bend; offset = o.offset_bend;

% number of vesicles
nv = numel(sortIdx(1,:));

Ypred = zeros(nv,o.nCompRelax);
% predict PCA coefficients for Xnew
Ypred(:,1:16) = predict(nets{1},Xinput(1:16,:,:,:));
% we use normalized output so take that back
Ypred(:,1:16) = (Ypred(:,1:16)-offset(1))*sdevChan1(1)/scale(1)+muChan1(1); 
if o.nCompRelax > 16
  % predict PCA coefficients for Xnew
  Ypred(:,17:32) = predict(nets{2},Xinput(17:32,:,:,:));
  % we use normalized output so take that back
  Ypred(:,17:32) = (Ypred(:,17:32)-offset(2))*sdevChan1(2)/scale(2)+muChan1(2); 
end

Xnew = zeros(2*N,nv);
for k = 1 : nv
  % reconstruct Xnew using PCA basis
  Xpred = (Ypred(k,:)*o.evects(:,1:o.nCompRelax)'+o.colMeans)';
    
  % destandardize
  Xpred = o.destandardize(Xpred,trans(:,k),rotate(k),scaling(k),sortIdx(:,k));

  % downsample to N
  Xnew(:,k) = [interpft(Xpred(1:end/2),N);interpft(Xpred(end/2+1:end),N)];
end
end % relaxWNN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xnew = translateVinfwNN(o,Xinput,vinf,Xold,rotate,sortIdx)
% Xinput is equally distributed in arc-length
% Xold as well. So, we add up coordinates of the same points.
N = numel(Xold(:,1))/2;
nv = numel(Xold(1,:));
Nnet = numel(sortIdx(:,1));

% load network files
FCnets = o.MVnets; activeModes = o.velActiveModes; outputSize = o.MVoutSize;

% Approximate the multiplication M*(FFTBasis)     
Z11r = zeros(Nnet,numel(activeModes),nv); Z12r = Z11r;
Z21r = Z11r; Z22r = Z11r;

for k = 1 : numel(activeModes)
  pred = predict(FCnets{k},Xinput)';
  Z11r(:,k,:) = interpft(pred(1:outputSize/4,:),Nnet);
  Z21r(:,k,:) = interpft(pred(outputSize/4+1:outputSize/2,:),Nnet);
  Z12r(:,k,:) = interpft(pred(outputSize/2+1:3*outputSize/4,:),Nnet);
  Z22r(:,k,:) = interpft(pred(3*outputSize/4+1:outputSize,:),Nnet);
end
% upsample vinf
vinfUp = [interpft(vinf(1:end/2,:),Nnet);interpft(vinf(end/2+1:end,:),Nnet)];
% Take fft of the velocity (should be standardized velocity)

MVinfMat = zeros(2*N,nv);
for k = 1 : nv
  % only sort points and rotate to pi/2 (no translation, no scaling)
  vinfStand = o.standardize(vinfUp(:,k),[0;0],rotate(k),1,sortIdx(:,k));
  z = vinfStand(1:end/2)+1i*vinfStand(end/2+1:end);

  zh = fft(z);
  V1 = real(zh(activeModes)); V2 = imag(zh(activeModes));
  % Compute the approximate value of the term M*vinf
  MVinfFull = [Z11r(:,:,k)*V1+Z12r(:,:,k)*V2; Z21r(:,:,k)*V1+Z22r(:,:,k)*V2];
  % Need to destandardize MVinf (take sorting and rotation back)
  MVinf = zeros(size(MVinfFull));
  MVinf([sortIdx(:,k);sortIdx(:,k)+Nnet]) = MVinfFull;
  MVinf = o.rotationOperator(MVinf,-rotate(k));
  % downsample MVinf
  MVinfMat(:,k) = [interpft(MVinf(1:end/2),N);interpft(MVinf(end/2+1:end),N)];
end
% Update the position
Xnew = Xold + o.dt*vinf-o.dt*MVinfMat;
end % translateVinfwNN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xnew = relaxWNNvariableKbDt(o,Xmid)
Nnet = o.Nnet;
N = numel(Xmid(:,1))/2;

% load network
nets = o.bendNets; muChan1 = o.muChan_bend; sdevChan1 = o.sdevChan_bend;
scale = o.scale_bend; offset = o.offset_bend;
KbDts = o.KbDts; flowKbDt = o.dt * o.kappa;

% number of vesicles
nv = numel(Xmid(1,:));

% number of nets used
nnets = numel(KbDts);

% Get 5 approximations, then interpolate then
scaling = zeros(nv,1); rotate = zeros(nv,1);
trans = zeros(2,nv); sortIdx = zeros(Nnet,nv);

for k = 1 : nv
  % 2) RELAXATION w/ NETWORK
  % Standardize vesicle Xmid
  [Xstand,scaling(k),rotate(k),trans(:,k),sortIdx(:,k)] = ...
      o.standardizationStep(Xmid(:,k));
  Xin = (Xstand'-o.colMeans)*o.evects(:,1:o.nCompRelax);
  
  for inet = 1 : nnets
    Xinput{inet}(1:16,1,1,k) = scale(inet,1)*(Xin(1:16)-muChan1(inet,1))/...
      sdevChan1(inet,1)+offset(inet,1);
    if o.nCompRelax > 16
    Xinput{inet}(17:32,1,1,k) = scale(inet,2)*(Xin(17:32)-muChan1(inet,2))/...
      sdevChan1(inet,2)+offset(inet,2);
    end
  end % inet = 1 : nnets
end % k = 1 : nv

Ypred = zeros(nv,o.nCompRelax,nnets);
for inet = 1 : nnets
% predict PCA coefficients for Xnew
Ypred(:,1:16,inet) = predict(nets{inet,1},Xinput{inet}(1:16,:,:,:));
% we use normalized output so take that back
Ypred(:,1:16,inet) = (Ypred(:,1:16,inet)-offset(inet,1))*sdevChan1(inet,1)/...
    scale(inet,1)+muChan1(inet,1); 
if o.nCompRelax > 16
  % predict PCA coefficients for Xnew
  Ypred(:,17:32,inet) = predict(nets{inet,2},Xinput{inet}(17:32,:,:,:));
  % we use normalized output so take that back
  Ypred(:,17:32,inet) = (Ypred(:,17:32,inet)-offset(inet,2))*sdevChan1(inet,2)/...
      scale(inet,2)+muChan1(inet,2); 
end
end

% Build Lagrange Interpolation Function
KbDts = log(KbDts); flowKbDt = log(flowKbDt);
funcl = zeros(nnets,1);
for inet = 1 : nnets
  xm = KbDts(inet); xj = [(1:inet-1)';(inet+1:nnets)'];
  funcl(inet) = prod(flowKbDt-KbDts(xj))/prod(xm-KbDts(xj));
end

Xnew = zeros(2*N,nv);
for k = 1 : nv
  YpredInt = Ypred(k,:,1)*funcl(1);
  for inet = 2 : nnets
    YpredInt = YpredInt + Ypred(k,:,inet)*funcl(inet);
  end
  % reconstruct Xnew using PCA basis
  Xpred = (YpredInt*o.evects(:,1:o.nCompRelax)'+o.colMeans)';
    
  % destandardize
  Xpred = o.destandardize(Xpred,trans(:,k),rotate(k),scaling(k),sortIdx(:,k));

  % downsample to N
  Xnew(:,k) = [interpft(Xpred(1:end/2),N);interpft(Xpred(end/2+1:end),N)];
end
end % relaxWNNvariableKbDt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xinput = prepareInputForNet(o,X,netType)
% Normalize input
if strcmp(netType,'tensionOnFourier')
scale = o.scale_tenV;
muChan = o.muChan_tenV;
sdevChan = o.sdevChan_tenV;
offset = o.offset_tenV;
elseif strcmp(netType,'tensionSelfBend')
scale = o.scale_tenBend;
muChan = o.muChan_tenBend;
sdevChan = o.sdevChan_tenBend;
offset = o.offset_tenBend;    
elseif strcmp(netType,'advection')
scale = o.scale_MV;
muChan = o.muChan_MV;
sdevChan = o.sdevChan_MV;
offset = o.offset_MV;    
elseif strcmp(netType,'relaxation')
scale = o.scale_bend;
muChan = o.muChan_bend;
sdevChan = o.sdevChan_bend;
offset = o.offset_bend;     
end

if strcmp(netType,'relaxation')
  % find PCA coefficients  
  Xinput(:,1,1,1) = (X'-o.colMeans)*o.evects(:,1:o.nCompRelax);  
  Xinput(1:16,1,1,1) = scale(1)*(Xinput(1:16,1,1,1)-muChan(1))/...
    sdevChan(1)+offset(1);
  if o.nCompRelax > 16
    Xinput(17:32,1,1,1) = scale(2)*(Xinput(17:32,1,1,1)-muChan(2))/...
        sdevChan(2)+offset(2);
  end
  if o.nCompRelax > 32
    Xinput(33:48,1,1,1) = scale(3)*(Xinput(33:48,1,1,1)-muChan(3))/...
        sdevChan(3)+offset(3);
  end
  if o.nCompRelax > 48
    Xinput(49:64,1,1,1) = scale(4)*(Xinput(49:64,1,1,1)-muChan(4))/...
        sdevChan(4)+offset(4);
  end
else
  % find PCA coefficients  
  Xinput(:,1,1,1) = (X'-o.colMeans)*o.evects(:,1:o.nComp);  
  Xinput(:,1,1,1) = scale*(Xinput(:,1,1,1)-muChan)/sdevChan+offset;        
end     
end % prepareInputForNet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,scaling,rotate,trans,sortIdx] = standardizationStep(o,Xin)

Nnet = o.Nnet;
oc = o.oc;

N = numel(Xin)/2;
if Nnet ~= N
  Xin = [interpft(Xin(1:end/2),Nnet);interpft(Xin(end/2+1:end),Nnet)];    
end

X = Xin;
% Equally distribute points in arc-length
for iter = 1 : 5
  [X,~,~] = oc.redistributeArcLength(X);
end
% Fix misalignment in center and angle due to reparametrization
X = oc.alignCenterAngle(Xin,X);

% standardize angle, center, scaling and point order
[trans,rotate,scaling,sortIdx] = o.referenceValues(X);
X = o.standardize(X,trans,rotate,scaling,sortIdx);
end % standardizationStep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function XrotSort = standardize(o,X,translation,rotation,scaling,sortIdx)
N = numel(sortIdx);

% translate, rotate and scale configuration
Xrotated = scaling*o.rotationOperator(o.translateOp(X,translation),rotation);   

% now order the points
XrotSort = [Xrotated(sortIdx);Xrotated(sortIdx+N)];

end % standardize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = destandardize(o,XrotSort,translation,rotation,scaling,sortIdx)

N = numel(sortIdx);    
    
% change ordering back 
X = zeros(size(XrotSort));
X([sortIdx;sortIdx+N]) = XrotSort;

% scaling back
X = X/scaling;

% take rotation back
cx = mean(X(1:end/2)); cy = mean(X(end/2+1:end));
X = o.rotationOperator([X(1:end/2)-cx;X(end/2+1:end)-cy],-rotation);
X = [X(1:end/2)+cx; X(end/2+1:end)+cy];

% take translation back
X = o.translateOp(X,-translation);

end % destandardize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [translation,rotation,scaling,sortIdx] = referenceValues(o,Xref)
oc = o.oc;
N = numel(Xref)/2;

% find translation, rotation and scaling
translation = [-mean(Xref(1:end/2));-mean(Xref(end/2+1:end))];
rotation = pi/2-oc.getIncAngle2(Xref);
    
% amount of scaling
[~,~,length] = oc.geomProp(Xref);
scaling = 1/length;
    
% find the ordering of the points
Xref = scaling*o.rotationOperator(o.translateOp(Xref,translation),rotation);

firstQuad = find(Xref(1:end/2)>=0 & Xref(end/2+1:end)>=0);
theta = atan2(Xref(end/2+1:end),Xref(1:end/2));
[~,idx]= min(theta(firstQuad));
sortIdx = [(firstQuad(idx):N)';(1:firstQuad(idx)-1)'];

end % referenceValues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xrot = rotationOperator(o,X,theta)
% Get x-y coordinates
Xrot = zeros(size(X));
x = X(1:end/2); y = X(end/2+1:end);

% Rotated shape
xrot = (x)*cos(theta) - (y)*sin(theta);
yrot = (x)*sin(theta) + (y)*cos(theta);

Xrot(1:end/2) = xrot;
Xrot(end/2+1:end) = yrot;
end % rotationOperator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xnew = translateOp(o,X,transXY)
Xnew = zeros(size(X));
Xnew(1:end/2) = X(1:end/2)+transXY(1);
Xnew(end/2+1:end) = X(end/2+1:end)+transXY(2);
end  % translateOp  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xnew,tenNew,etaNew,etaIntNew,etaExtNew,RSnew] = DNNsolveMixture(o,...
        Xold,tenOld,etaOld,etaIntOld,etaExtOld,RSold)

% this is used to take a time step with mixtures of particles (e.g., 
% vesicles of a certain reduced area and circular particles). Circular
% particles should be the first 'nCircular' columns of Xold, Xnew. 
    
tt = o.tt;
op = tt.op;
if o.confined
  if isempty(o.wallsInt) % not DLD
  opWall = tt.opWall;
  walls = o.walls;
  uWalls = walls.u;
  nvbd = walls.nv;
  Nbd = walls.N;
  diffDiscWalls = false;
  etaIntNew = []; etaExtNew = [];
  else
  opWallInt = tt.opWallInt;
  opWallExt = tt.opWallExt;  
  wallsInt = o.wallsInt;
  wallsExt = o.wallsExt;
  nvbdInt = wallsInt.nv;
  nvbdExt = 1;
  nvbd = nvbdExt + nvbdInt;
  NbdInt = wallsInt.N;
  NbdExt = wallsExt.N;
  diffDiscWalls = true;
  etaNew = [];
  end % if isempty(o.wallsInt)
  RSnew = zeros(3,nvbd);
else
  vback = o.vinf(Xold);
  etaNew = []; RSnew = []; etaIntNew = []; etaExtNew = [];
end

% explicit time stepping w/ splitting
vesicle = capsules(Xold,[],[],o.kappa,1,1);
vesicle.setUpRate();
nv = vesicle.nv; % this includes nCircular
N = vesicle.N;
nCircular = o.nCircular;

if tt.fmm
  Nup = op.LPuprate * vesicle.N;
  Xup = [interpft(Xold(1:end/2,:),Nup);interpft(Xold(end/2+1:end,:),Nup)];
  vesicleUp = capsules(Xup,[],[],vesicle.kappa,vesicle.viscCont,0); 
  SLPnoCorr = op.stokesSLmatrixNoCorr(vesicleUp);
else
  SLPnoCorr = [];
end

% Get the near structure
if o.confined
  if ~diffDiscWalls  
    [NearV2V,NearV2W] = vesicle.getZone(walls,3);
    [~,NearW2V] = walls.getZone(vesicle,2);
    if nvbd > 1
      NearW2W = o.NearW2W;  
    end % if nvbd > 1
  else
    % Get the near structure, so that we know what to neglect
    [NearV2V,NearV2Wint] = vesicle.getZone(wallsInt,3);
      
    if isempty(tt.NearWint2Wint)
      [tt.NearWint2Wint,NearWint2V] = wallsInt.getZone(vesicle,3);
    else
      % there is no need to compute W2W again, since they do not move
      [~,NearWint2V] = wallsInt.getZone(vesicle,2);
    end
    NearWint2Wint = tt.NearWint2Wint;

    [~,NearV2Wext] = vesicle.getZone(wallsExt,2);
    [~,NearWext2V] = wallsExt.getZone(vesicle,2);

    if isempty(tt.NearWint2Wext)
      % there is no need to compute W2W again, since they do not move
      [~,tt.NearWint2Wext] = wallsInt.getZone(wallsExt,2);
    end
    NearWint2Wext = tt.NearWint2Wext;

    if isempty(tt.NearWext2Wint)
      % there is no need to compute W2W again, since they do not move
      [~,tt.NearWext2Wint] = wallsExt.getZone(wallsInt,2);    
    end
    NearWext2Wint = tt.NearWext2Wint;
  end % if ~diffDiscWalls
else
  NearV2V = vesicle.getZone([],1);    
end % if o.confined


% 1) EXPLICIT TENSION AT THE CURRENT STEP

% Compute bending forces + old tension forces
fBend = vesicle.tracJump(Xold,zeros(N,nv));
fTen = vesicle.tracJump(zeros(2*N,nv),tenOld);
tracJump = fBend+fTen;

if nv > 1
  % SLP due to traction jump
  kernelRegul = @op.exactStokesSLregul; % regularized kernel for near-field
  kernelDirect = @op.exactStokesSL; % direct evaluation for far-field
  if tt.fmm
    kernel = @op.exactStokesSLfmm;
  else
    kernel = @op.exactStokesSLregul;
  end
  ves2vesSLP = op.evaluateLayerPoten(vesicle,tracJump,SLPnoCorr,NearV2V,...
        kernel,kernelDirect,kernelRegul,vesicle,true);
else
  ves2vesSLP = zeros(2*N,1);
end

% DLP due to walls or background flow
if o.confined
  if ~diffDiscWalls  
    kernelRegul = @opWall.exactStokesDLregul;
    kernelDirect = @opWall.exactStokesDL;
    if tt.fmmDLP
      kernel = @opWall.exactStokesDLnewfmm;
    else
      kernel = @opWall.exactStokesDLregul;
    end
    vback = opWall.evaluateLayerPoten(walls,etaOld,[],NearW2V,...
        kernel,kernelDirect,kernelRegul,vesicle,false);  
  else % then DLD
    kernelRegul = @opWallInt.exactStokesDLregul;
    kernelRegul2 = @opWallExt.exactStokesDLregul;
    kernelDirect = @opWallInt.exactStokesDL;
    kernelDirect2 = @opWallExt.exactStokesDL;
    if tt.fmmDLP
      kernel = @opWallInt.exactStokesDLnewfmm;
      kernel2 = @opWallExt.exactStokesDLnewfmm;
    else
      kernel = @opWallInt.exactStokesDLregul;
      kernel2 = @opWallExt.exactStokesDLregul;
    end
    % interior walls to vesicles
    wallInt2vesDLP = opWallInt.evaluateLayerPoten(wallsInt,etaIntOld,[],NearWint2V,...
        kernel,kernelDirect,kernelRegul,vesicle,false);
    % exterior wall to vesicles
    wallExt2vesDLP = opWallExt.evaluateLayerPoten(wallsExt,etaExtOld,[],NearWext2V,...
        kernel2,kernelDirect2,kernelRegul2,vesicle,false);
    
    vback = wallInt2vesDLP + wallExt2vesDLP;
  end %if ~diffDiscWalls
  
  for k = 2:nvbd
    stokeslet = RSold(1:2,k); rotlet = RSold(3,k);  
    if ~diffDiscWalls
      vback = vback + tt.RSlets(vesicle.X,walls.center(:,k),stokeslet,rotlet);
    else
      vback = vback + tt.RSlets(vesicle.X,wallsInt.center(:,k-1),stokeslet,rotlet);
    end
  end % k = 2:nvbd (stokeslet and rotlet)
  
end % if o.confined

scaling = zeros(nv-nCircular,1); rotate = zeros(nv-nCircular,1); 
trans = zeros(2,nv-nCircular);
sortIdx = zeros(o.Nnet,nv-nCircular);
XinputFourier = zeros(o.nComp,1,1,nv-nCircular);
XinputTenBen = zeros(o.nComp,1,1,nv-nCircular);
for k = nCircular+1 : nv
  % Standardize vesicle  
  [Xstand,scaling(k-nCircular),rotate(k-nCircular),trans(:,k-nCircular),...
      sortIdx(:,k-nCircular)] = o.standardizationStep(Xold(:,k));
  % Prepare inputs for Fourier network and self-bending network
  XinputFourier(:,:,:,k-nCircular) = o.prepareInputForNet(Xstand,'tensionOnFourier');
  XinputTenBen(:,:,:,k-nCircular) = o.prepareInputForNet(Xstand,'tensionSelfBend');
end
% Approximate inv(Div*G*Ten)*Div*vExt
vBackSolve = o.invTenMatOnVback(XinputFourier,vback(:,nCircular+1:end)...
    +ves2vesSLP(:,nCircular+1:end),scaling,rotate,sortIdx);
% Approximate inv(Div*G*Ten)*G*(-Ben)*x
selfBendSolve = o.invTenMatOnSelfBend(XinputTenBen,scaling,sortIdx,N);
% update tension
tenNew = -(vBackSolve+selfBendSolve);


% Evaluate tension for circular particles
tenCirc = zeros(N,nCircular); 
for k = 1 : nCircular
  [XcircStand,scalingCirc,rotateCirc,~,sortIdxCirc] = ...
      o.standardizationStep(Xold(:,k));
  vOnCirc = vback(:,k)+ves2vesSLP(:,k);
  vOnCirc = [interpft(vOnCirc(1:end/2),o.Nnet);interpft(vOnCirc(end/2+1:end),o.Nnet)];
  vOnCircStand = o.standardize(vOnCirc,[0;0],rotateCirc,1,sortIdxCirc);
  tenOthersStand = o.invTenMatCirc*vOnCircStand;
  tenOthers = zeros(size(tenOthersStand));
  tenOthers(sortIdxCirc) = tenOthersStand;
  
  tenSelfStand = o.selfBendMatCirc*XcircStand;
  tenSelf = zeros(size(tenSelfStand));
  tenSelf(sortIdxCirc) = tenSelfStand/scalingCirc^2;
  tenCirc(:,k) = interpft(-tenOthers-tenSelf,N);
end
tenNew = [tenCirc tenNew];
  
% Update the traction jump
fTen = vesicle.tracJump(zeros(2*N,nv),tenNew); tracJump = fBend+fTen;

% 2) SOLVE FOR DENSITY and RS ON WALLS and UPDATE wall2vesDLP
if o.confined
  % vesicle2wall interactions
  kernelRegul = @op.exactStokesSLregul;
  kernelDirect = @op.exactStokesSL;
  if tt.fmm
    kernel = @op.exactStokesSLfmm;
  else
    kernel = @op.exactStokesSLregul;
  end
  
  % SLP due to vesicles on walls
  if ~diffDiscWalls
    ves2wallSLP = op.evaluateLayerPoten(vesicle,tracJump,[],NearV2W,...
        kernel,kernelDirect,kernelRegul,walls,0);
  
    RHS = [uWalls(:)-ves2wallSLP(:);zeros(3*(nvbd-1),1)];
    etaRS = tt.bdiagWall*RHS;
    etaNew = zeros(2*Nbd,nvbd);
    for iw = 1 : nvbd
      etaNew(:,iw) = etaRS((iw-1)*2*Nbd+1:iw*2*Nbd);  
      if iw <= nvbd-1
        RSnew(:,iw+1) = etaRS(2*Nbd*nvbd+(iw-1)*3+1:2*Nbd*nvbd+iw*3);
      end
    end 
  else % if DLD  
    ves2wallInt = op.evaluateLayerPoten(vesicle,tracJump,[],NearV2Wint,...
        kernel,kernelDirect,kernelRegul,wallsInt,0);
    ves2wallExt = op.evaluateLayerPoten(vesicle,tracJump,[],NearV2Wext,...
        kernel,kernelDirect,kernelRegul,wallsExt,0);
    
    RHSint = wallsInt.u(:)-ves2wallInt(:);
    RHSext = wallsExt.u(:)-ves2wallExt(:);
    RHS = [RHSext;RHSint;zeros(3*(nvbd-1),1)];
    etaRS = tt.bdiagWall*RHS;
    etaExtNew = etaRS(1:2*NbdExt);
    etaRS = etaRS(2*NbdExt+1:end);
    etaIntNew = zeros(2*NbdInt,nvbdInt);
    for iw = 1 : nvbdInt
      etaIntNew(:,iw) = etaRS((iw-1)*2*NbdInt+1:iw*2*NbdInt);  
      if iw <= nvbdInt
        RSnew(:,iw+1) = etaRS(2*NbdInt*nvbdInt+(iw-1)*3+1:2*NbdInt*nvbdInt+iw*3);
      end
    end   
  end % if ~diffDiscWalls
  
  % Update vback due to new density,rotlets and stokeslet
  if ~diffDiscWalls  
    kernelRegul = @opWall.exactStokesDLregul;
    kernelDirect = @opWall.exactStokesDL;
    if tt.fmmDLP
      kernel = @opWall.exactStokesDLnewfmm;
    else
      kernel = @opWall.exactStokesDLregul;
    end
    vback = opWall.evaluateLayerPoten(walls,etaNew,[],NearW2V,...
        kernel,kernelDirect,kernelRegul,vesicle,false);  
  else % then DLD
    kernelRegul = @opWallInt.exactStokesDLregul;
    kernelRegul2 = @opWallExt.exactStokesDLregul;
    kernelDirect = @opWallInt.exactStokesDL;
    kernelDirect2 = @opWallExt.exactStokesDL;
    if tt.fmmDLP
      kernel = @opWallInt.exactStokesDLnewfmm;
      kernel2 = @opWallExt.exactStokesDLnewfmm;
    else
      kernel = @opWallInt.exactStokesDLregul;
      kernel2 = @opWallExt.exactStokesDLregul;
    end
    % interior walls to vesicles
    wallInt2vesDLP = opWallInt.evaluateLayerPoten(wallsInt,etaIntNew,[],NearWint2V,...
        kernel,kernelDirect,kernelRegul,vesicle,false);
    % exterior wall to vesicles
    wallExt2vesDLP = opWallExt.evaluateLayerPoten(wallsExt,etaExtNew,[],NearWext2V,...
        kernel2,kernelDirect2,kernelRegul2,vesicle,false);
    
    vback = wallInt2vesDLP + wallExt2vesDLP;
  end %if ~diffDiscWalls
  
  for k = 2:nvbd
    stokeslet = RSnew(1:2,k); rotlet = RSnew(3,k);  
    if ~diffDiscWalls
      vback = vback + tt.RSlets(vesicle.X,walls.center(:,k),stokeslet,rotlet);
    else
      vback = vback + tt.RSlets(vesicle.X,wallsInt.center(:,k-1),stokeslet,rotlet);
    end
  end % k = 2:nvbd (stokeslet and rotlet)
end % if iconfined    

% 3) SOLVE FOR NEXT TIME STEP USING OPERATOR SPLITTING
% EXPLICIT-TREATMENT OF INTER-VESICLE BENDING AND TENSION
% TENSION AND X are coupled (semi-implicit).
if nv > 1
  % SLP due to traction jump
  kernelRegul = @op.exactStokesSLregul; % regularized kernel for near-field
  kernelDirect = @op.exactStokesSL; % direct evaluation for far-field
  if tt.fmm
    kernel = @op.exactStokesSLfmm;
  else
    kernel = @op.exactStokesSLregul;
  end
  ves2vesSLP = op.evaluateLayerPoten(vesicle,tracJump,SLPnoCorr,NearV2V,...
        kernel,kernelDirect,kernelRegul,vesicle,true);
else
  ves2vesSLP = zeros(2*N,1);
end

% Walls2Ves is already computed if confined (vback)    
rotate = zeros(nv-nCircular,1); trans = zeros(2,nv-nCircular); 
sortIdx = zeros(o.Nnet,nv-nCircular);
scaling = zeros(nv-nCircular,1); Xinput = zeros(o.nComp,1,1,nv-nCircular);
for k = nCircular+1 : nv
  % 1) TRANSLATION W/ NETWORK
  % Standardize vesicle
  [Xstand,scaling(k-nCircular),rotate(k-nCircular),trans(:,k-nCircular),...
      sortIdx(:,k-nCircular)] = o.standardizationStep(Xold(:,k));
  % Prepare input for advection network
  Xinput(:,1,1,k-nCircular) = o.prepareInputForNet(Xstand,'advection');
end
% Take a step due to advaction
Xmid = o.translateVinfwNN(Xinput,vback(:,nCircular+1:end) ...
    + ves2vesSLP(:,nCircular+1:end),Xold(:,nCircular+1:end),rotate,sortIdx);

% Advection step for circular particles
XmidCirc = zeros(2*N,nCircular);
for k = 1 : nCircular
  [~,~,rotateCirc,~,sortIdxCirc] = o.standardizationStep(Xold(:,k));
  vOnCirc = vback(:,k)+ves2vesSLP(:,k);
  vOnCircUp = [interpft(vOnCirc(1:end/2),o.Nnet);...
      interpft(vOnCirc(end/2+1:end),o.Nnet)];
  vOnCircStand = o.standardize(vOnCircUp,[0;0],rotateCirc,1,sortIdxCirc);
  MVinfCircStand = o.Mcirc*vOnCircStand;
  MVinfCirc = zeros(size(MVinfCircStand));
  MVinfCirc([sortIdxCirc;sortIdxCirc+o.Nnet]) = MVinfCircStand;
  MVinfCirc = o.rotationOperator(MVinfCirc,-rotateCirc);
  MVinfCirc = [interpft(MVinfCirc(1:end/2),N);interpft(MVinfCirc(end/2+1:end),N)];
  XmidCirc(:,k) = Xold(:,k) + o.dt*vOnCirc - o.dt*MVinfCirc;
end
Xmid = [XmidCirc Xmid];

if ~o.variableKbDt
  rotate = zeros(nv-nCircular,1); trans = zeros(2,nv-nCircular); 
  sortIdx = zeros(o.Nnet,nv-nCircular); scaling = zeros(nv-nCircular,1); 
  Xinput = zeros(o.nCompRelax,1,1,nv-nCircular);
  for k = nCircular+1 : nv
    % 2) RELAXATION w/ NETWORK
    % Standardize vesicle Xmid
    [Xstand,scaling(k-nCircular),rotate(k-nCircular),trans(:,k-nCircular),...
        sortIdx(:,k-nCircular)] = o.standardizationStep(Xmid(:,k));
    % Prepare input for relaxation network
    Xinput(:,1,1,k-nCircular) = o.prepareInputForNet(Xstand,'relaxation');
  end
  % Take a step
  Xnew = o.relaxWNN(Xinput,scaling,rotate,trans,sortIdx,N);
else
  Xnew = o.relaxWNNvariableKbDt(Xmid(:,nCircular+1:end));    
end

% Relaxation for circular particles
XnewCirc = zeros(2*N,nCircular);
for k = 1 : nCircular
  [Xstand,scalingCirc,rotateCirc,transCirc,sortIdxCirc] = ...
      o.standardizationStep(Xmid(:,k));
  XnewStand = o.relaxMatCirc*Xstand;
  XnewUp = o.destandardize(XnewStand,transCirc,rotateCirc,scalingCirc,sortIdxCirc);
  XnewCirc(:,k) = [interpft(XnewUp(1:end/2),N);interpft(XnewUp(end/2+1:end),N)];
end
Xnew = [XnewCirc Xnew];

end % DNNsolveMixture
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xnew = relaxExactSolve(o,vesicle,vinf,dt,Xold,op)
% HERE, BENDING IN TENSION SOLVE IS IMPLICIT

% SLP
G = op.stokesSLmatrix(vesicle);
% Bending, tension and surface divergence
[Ben,Ten,Div] = vesicle.computeDerivs;
M = G*Ten*((Div*G*Ten)\eye(vesicle.N))*Div;
rhs = Xold + dt*(eye(2*vesicle.N)-M)*vinf;
LHS = (eye(2*vesicle.N)-vesicle.kappa*dt*(-G*Ben+M*G*Ben));
Xnew = LHS\rhs;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xnew = translateVinf(o,vinf,dt,Xold,op)
% ADVECTION PART FOR OPERATOR SPLITTING
vesicle = capsules(Xold,[],[],1,1,1);
vesicle.setUpRate();
 
G = op.stokesSLmatrix(vesicle);
[~,Ten,Div] = vesicle.computeDerivs;

M = G*Ten*((Div*G*Ten)\eye(vesicle.N))*Div;
Xnew = Xold + dt*(eye(2*vesicle.N)-M)*vinf;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vinf = setBgFlow(o,bgFlow,speed)

vinf = @(X) zeros(size(X));      
if strcmp(bgFlow,'relax')
  vinf = @(X) zeros(size(X));  % relaxation
elseif strcmp(bgFlow,'shear') 
  vinf = @(X) speed*[X(end/2+1:end,:);zeros(size(X(1:end/2,:)))]; 
elseif strcmp(bgFlow,'tayGreen')
  vinf = @(X) speed*[sin(X(1:end/2,:)).*cos(X(end/2+1:end,:));-...
    cos(X(1:end/2,:)).*sin(X(end/2+1:end,:))]; % Taylor-Green
elseif strcmp(bgFlow,'parabolic')
  vinf = @(X) [speed*(1-(X(end/2+1:end,:)/0.2).^2);...
      zeros(size(X(1:end/2,:)))];
elseif strcmp(bgFlow,'rotation')
  vinf = @(X) [-sin(atan2(X(end/2+1:end,:),X(1:end/2,:)))./sqrt(X(1:end/2,:).^2+X(end/2+1:end,:).^2);...
    cos(atan2(X(end/2+1:end,:),X(1:end/2,:)))./sqrt(X(1:end/2,:).^2+X(end/2+1:end,:).^2)]*speed;
elseif strcmp(bgFlow,'extensional')
  vinf = @(X) speed*[-X(1:end/2,:);X(end/2+1:end,:)];
end
    
end % setBgFlow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tt = buildTstep(o,X,pramsIn)
% Assign tstep class' parameters
prams.N = pramsIn.N;
prams.nv = pramsIn.nv;
prams.totnv = pramsIn.totnv;
prams.Nbd = pramsIn.Nbd;

prams.NbdInt = pramsIn.NbdInt;
prams.NbdExt = pramsIn.NbdExt;
prams.nvbd = max(pramsIn.nvbdInt+pramsIn.nvbdExt,pramsIn.nvbd);
prams.nvbdInt = pramsIn.nvbdInt; 
prams.nvbdExt = pramsIn.nvbdExt;

prams.nrow = pramsIn.nrow;
prams.ncol = pramsIn.ncol;
prams.Dpostx = pramsIn.Dpostx;
prams.Dposty = pramsIn.Dposty;
prams.Dx = pramsIn.Dx; 
prams.Dy = pramsIn.Dy;
prams.epsilon = pramsIn.epsilon;

options.diffDiscWalls = 0;
if prams.nvbdInt > 0
  options.diffDiscWalls = 1;
end

prams.kappa = pramsIn.kappa;

options.verbose = 0;
options.saveData = 0;
options.usePlot = 0;
options.track = 0;
options.quiver = 0;
options.axis = 0;
prams.T = pramsIn.Th;
prams.gmresTol = 1e-8;
prams.m = pramsIn.Th/pramsIn.dt;
prams.errTol = 1e-1;
options.tracers = 0;
options.timeAdap = 0;
options.order = 1;
prams.areaLenTol = 1e-1;
options.fmm = pramsIn.fmm;
options.fmmDLP = pramsIn.fmmDLP;

options.correctShape = true;
options.adhesion = 0;
options.repulsion = 0;

options.confined = o.confined;
options.farField = pramsIn.bgFlow;
options.farFieldSpeed = pramsIn.speed;

[options,prams] = initVes2D(options,prams);
    
om = monitor(X,options,prams);
tt = tstep(options,prams,om);
if ~options.confined % if not confined, then assign vinf
  tt.farField = o.setBgFlow(pramsIn.bgFlow,pramsIn.speed); 
end

if strcmp(pramsIn.bgFlow,'rotDLD')
  tt.farField = @(X,Xint) tt.bgFlow(X,options.farField,'intWalls',Xint,...
      'velG',pramsIn.gPer,'Speed',1);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % methods

end % dnnTools
