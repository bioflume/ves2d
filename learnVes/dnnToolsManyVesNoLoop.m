classdef dnnToolsManyVesNoLoop

properties
KbDts   
repLenScale
variableKbDt
confined
ttExact
tt
dt
vinf
walls
NearW2W
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
oc
useNear
useTrueNear
repulsion
repStrength
minDist
kappa
interpOrder
Nexact
N
dtExact
end

methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = dnnToolsManyVesNoLoop(X,Xwalls,prams)
if nargin > 0    
if strcmp(prams.bgFlow,'couette')    
o.confined = true;
else
o.confined = false;
end
o.interpOrder = prams.interpOrder;
o.kappa = prams.kappa;
o.dt = prams.dt;    
o.dtExact = prams.dtExact;
o.Nexact = prams.Nexact;
o.N = prams.N;
prams.N = o.N;
prams.dt = o.dt;
o.tt = o.buildTstep(X,prams); 
prams.N = o.Nexact;
prams.dt = o.dtExact;
o.ttExact = o.buildTstep(X,prams);
  if o.confined
    o = o.initWalls(Xwalls);
    o.vinf = [];
  else
    o.vinf = o.setBgFlow(prams.bgFlow,prams.speed);  
  end % if confined
end % if nargin > 0

end % dnnTools2Ves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = initWalls(o,Xwalls)
nvbd = numel(Xwalls(1,:)); % number of walls

tt = o.tt;
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
ttExact = o.ttExact;
ttExact.wallDLP = tt.wallDLP;
ttExact.wallDLPnoCorr = tt.wallDLPnoCorr;
% N0 to remove rank-1 deficiency
tt.wallN0 = opWall.stokesN0matrix(walls);
ttExact.wallN0 = tt.wallN0;

% inverse of the walls matrix (-1/2*density + DLP + N0 + R-S)
tt.matFreeWalls = true; % so, create walls matrix and save it
tt.bdiagWall = tt.wallsPrecond(walls); % inverse of the walls matrix
ttExact.matFreeWalls = tt.matFreeWalls;
ttExact.bdiagWall = tt.bdiagWall;
o.tt = tt;
o.ttExact = ttExact;
% tt includes tt.wallDLPandRSmat which is walls matrix saved 

end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xnew,tenNew,etaNew,RSnew,vbackTot] = DNNsolve(o,Xold,tenOld,etaOld,RSold,Nnet)
tt = o.tt;
op = tt.op;
if o.confined
  opWall = tt.opWall;
  invWallsMat = tt.bdiagWall;
  walls = o.walls;
  uWalls = walls.u;
  nvbd = walls.nv;
  Nbd = walls.N;
  etaNew = zeros(2*Nbd,nvbd);
  RSnew = zeros(3,nvbd);
else
  vback = o.vinf(Xold);
  etaNew = []; RSnew = [];
end

% semi-implicit time stepping w/ splitting
vesicle = capsules(Xold,[],[],o.kappa,1,1);
vesicle.setUpRate();
nv = vesicle.nv;
N = vesicle.N;
SLPdiag = op.stokesSLmatrixNoCorr(vesicle);
if tt.fmm
  Nup = op.LPuprate * vesicle.N;
  Xup = [interpft(Xold(1:end/2,:),Nup);interpft(Xold(end/2+1:end,:),Nup)];
  vesicleUp = capsules(Xup,[],[],vesicle.kappa,vesicle.viscCont,0); 
  SLPnoCorr = op.stokesSLmatrixNoCorr(vesicleUp);
else
  SLPnoCorr = [];
end

% Get the near structure, so that we know what to neglect
if o.confined
  [NearV2V,NearV2W] = vesicle.getZone(walls,3);
  [~,NearW2V] = walls.getZone(vesicle,2);
  if nvbd > 1
    NearW2W = o.NearW2W;  
  end % if nvbd > 1
else
  NearV2V = vesicle.getZone([],1);    
end % if o.confined

% Compute bending forces + old tension forces
fBend = vesicle.tracJump(Xold,zeros(N,nv));
fTen = vesicle.tracJump(zeros(2*N,nv),tenOld);
tracJump = fBend+fTen;

% 1) EXPLICIT TENSION AT THE CURRENT STEP
% neglect near interactions

if nv > 1
  % Far-field due to traction jump
  kernelDirect = @op.exactStokesSL;
  if tt.fmm
    kernel = @op.exactStokesSLfmm;
  else
    kernel = @op.exactStokesSL;
  end
  
  if o.useTrueNear
    tt.Galpert = op.stokesSLmatrix(vesicle);  
    SLP = @(X) op.exactStokesSLdiag(vesicle,tt.Galpert,X);
    farFieldtracJump = op.nearSingInt(vesicle,tracJump,SLP,SLPnoCorr,...
        NearV2V,kernel,kernelDirect,vesicle,true,false);
    % if repulsion is on, compute velocity due to that
    if o.repulsion
      if o.confined
        repulsion = vesicle.repulsionSchemeSimple(Xold,o.repStrength,o.repLenScale,...
            walls,[],[]);
      else
        repulsion = vesicle.repulsionSchemeSimple(Xold,o.repStrength,o.repLenScale,...
            [],[],[]);
      end
      Frepulsion = op.exactStokesSLdiag(vesicle,tt.Galpert,repulsion) + ...
          op.nearSingInt(vesicle,repulsion,SLP,SLPnoCorr,...
          NearV2V,kernel,kernelDirect,vesicle,true,false);
      farFieldtracJump = farFieldtracJump + Frepulsion; 
    end % if o.repulsion

  else
    [~,nearField,farFieldtracJump] = op.divideNearFarSLP(vesicle,...
        tracJump,SLPnoCorr,NearV2V,kernel,kernelDirect,vesicle,true); 
    if o.repulsion
      if o.confined
        repulsion = vesicle.repulsionSchemeSimple(Xold,o.repStrength,o.repLenScale,...
            walls,[],[]);
      else
        repulsion = vesicle.repulsionSchemeSimple(Xold,o.repStrength,o.repLenScale,...
            [],[],[]);
      end
      FrepSelf = zeros(2*vesicle.N, vesicle.nv);
      for iv = 1 : vesicle.nv
        FrepSelf(:,iv) = SLPdiag(:,:,iv) * repulsion(:,iv);    
      end
      [~,FrepNear,FrepFar] = op.divideNearFarSLP(vesicle, repulsion, SLPnoCorr, NearV2V, kernel, kernelDirect, vesicle,true);
      Frepulsion = FrepSelf + FrepNear + FrepFar;
    end
    if o.useNear
      farFieldtracJump = farFieldtracJump + Frepulsion + 0*nearField;
    end
  end
  
else
  farFieldtracJump = zeros(2*N,1);
end

% Far-field due to walls or background flow
if o.confined
  kernelDirect = @opWall.exactStokesDL;
  if tt.fmmDLP
    kernel = @opWall.exactStokesDLnewfmm;
  else
    kernel = @opWall.exactStokesDL;
  end
  
  if o.useTrueNear
    DLP = @(X) -1/2*X + opWall.exactStokesDLdiag(walls,tt.wallDLP,X); 
    vback = opWall.nearSingInt(walls,etaOld,DLP,[],NearW2V,kernel,kernelDirect,...
        vesicle,false,false);
  else
    [~,nearField,vback] = opWall.divideNearFarSLP(walls,etaOld,[],NearW2V,...
        kernel,kernelDirect,vesicle,false);
    if o.useNear
      vback = vback + 0*nearField;
    end
  end
  
  for k = 2:nvbd
    stokeslet = RSold(1:2,k); rotlet = RSold(3,k);  
    vback = vback + tt.RSlets(vesicle.X,walls.center(:,k),stokeslet,rotlet);
  end
end

scaling = zeros(nv,1); rotate = zeros(nv,1); trans = zeros(2,nv);
sortIdx = zeros(Nnet,nv);
XinputFourier = zeros(o.nComp,1,1,nv);
XinputTenBen = zeros(o.nComp,1,1,nv);
for k = 1 : nv
  % Standardize vesicle  
  [Xstand,scaling(k),rotate(k),trans(:,k),sortIdx(:,k)] = ...
      o.standardizationStep(Xold(:,k),Nnet);
  % Prepare inputs for Fourier network and self-bending network
  XinputFourier(:,:,:,k) = o.prepareInputForNet(Xstand,'tensionOnFourier');
  XinputTenBen(:,:,:,k) = o.prepareInputForNet(Xstand,'tensionSelfBend');
end
% Approximate inv(Div*G*Ten)*Div*vExt
vBackSolve = o.invTenMatOnVback(XinputFourier,vback+farFieldtracJump,...
    scaling,rotate,sortIdx);
% Approximate inv(Div*G*Ten)*G*(-Ben)*x
selfBendSolve = o.invTenMatOnSelfBend(XinputTenBen,scaling,sortIdx,N);
% update tension
tenNew = -(vBackSolve+selfBendSolve);

% Update the traction jump
fTen = vesicle.tracJump(zeros(2*N,nv),tenNew); tracJump = fBend+fTen;

if o.confined
  % 2) SOLVE FOR DENSITY and RS ON WALLS
  % vesicle2wall interactions
  kernelDirect = @op.exactStokesSL;
  if tt.fmm
    kernel = @op.exactStokesSLfmm;
  else
    kernel = @op.exactStokesSL;
  end
  
  if o.useTrueNear
    SLP = @(X) op.exactStokesSLdiag(vesicle,tt.Galpert,X);  
    ves2wallInt = op.nearSingInt(vesicle,tracJump,SLP,[],NearV2W,...
        kernel,kernelDirect,walls,false,false);
    if o.repulsion
      FrepWall = op.nearSingInt(vesicle,repulsion,SLP,[],NearV2W,kernel,...
          kernelDirect,walls,false,false);
      ves2wallInt = ves2wallInt + FrepWall;
    end % if o.repulsion
  else
    [~,nearField,ves2wallInt] = op.divideNearFarSLP(vesicle,tracJump,[],NearV2W,...
        kernel,kernelDirect,walls,0);
    if o.repulsion
      [~,FrepWallNear,FrepWallFar] = op.divideNearFarSLP(vesicle, repulsion, [], NearV2W, kernel, kernelDirect, walls,0);
      FrepulsionWall = FrepWallNear + FrepWallFar;
    end
    if o.useNear
      ves2wallInt = FrepulsionWall + ves2wallInt + 0*nearField;
    end
  end
  
  RHS = [uWalls(:)-ves2wallInt(:);zeros(3*(nvbd-1),1)];
  etaRS = invWallsMat*RHS;
  for iw = 1 : nvbd
    etaNew(:,iw) = etaRS((iw-1)*2*Nbd+1:iw*2*Nbd);  
    if iw <= nvbd-1
      RSnew(:,iw+1) = etaRS(2*Nbd*nvbd+(iw-1)*3+1:2*Nbd*nvbd+iw*3);
    end
  end 
  % Update vback due to new density,rotlets and stokeslet
  kernelDirect = @opWall.exactStokesDL;
  if tt.fmmDLP
    kernel = @opWall.exactStokesDLnewfmm;
  else
    kernel = @opWall.exactStokesDL;
  end
  
  if o.useTrueNear
    DLP = @(X) -1/2*X + opWall.exactStokesDLdiag(walls,tt.wallDLP,X);   
    vback = opWall.nearSingInt(walls,etaNew,DLP,[],NearW2V,kernel,...
        kernelDirect,vesicle,false,false);
  else
    [~,nearField,vback] = opWall.divideNearFarSLP(walls,etaNew,[],NearW2V,...
        kernel,kernelDirect,vesicle,false);
    if o.useNear
      vback = vback + 0*nearField;
    end
  end
  
  for k = 2:nvbd
    stokeslet = RSnew(1:2,k); rotlet = RSnew(3,k);  
    vback = vback + tt.RSlets(vesicle.X,walls.center(:,k),stokeslet,rotlet);
  end
end % if iconfined    

% 3) SOLVE FOR NEXT TIME STEP USING OPERATOR SPLITTING
% EXPLICIT-TREATMENT OF INTER-VESICLE BENDING AND TENSION
% TENSION AND X are coupled (semi-implicit).
Xnew = zeros(size(Xold));
if nv > 1
  % compute traction jump due to explicit tension and old bending
  kernelDirect = @op.exactStokesSL;
  if tt.fmm
    kernel = @op.exactStokesSLfmm;
  else
    kernel = @op.exactStokesSL;
  end
  
  if o.useTrueNear
    SLP = @(X) op.exactStokesSLdiag(vesicle,tt.Galpert,X);    
    farFieldtracJump = op.nearSingInt(vesicle,tracJump,SLP,SLPnoCorr,NearV2V,kernel,...
        kernelDirect,vesicle,true,false);
    % if repulsion is on, add the velocity (alread computed) due to that
    if o.repulsion
      farFieldtracJump = farFieldtracJump + Frepulsion; 
    end % if o.repulsion
  else
    [~,nearField,farFieldtracJump] = op.divideNearFarSLP(...
            vesicle,tracJump,SLPnoCorr,NearV2V,kernel,kernelDirect,vesicle,1);
    if o.useNear
      farFieldtracJump = farFieldtracJump + 0*nearField;
      if o.repulsion; farFieldtracJump = farFieldtracJump + Frepulsion; end;
    end    
  end
  
else
  farFieldtracJump = zeros(2*N,1);  
end

% Walls2Ves is already computed if confined (vback)    
rotate = zeros(nv,1); trans = zeros(2,nv); sortIdx = zeros(Nnet,nv);
scaling = zeros(nv,1); Xinput = zeros(o.nComp,1,1,nv);
for k = 1 : nv
  % 1) TRANSLATION W/ NETWORK
  % Standardize vesicle
  [Xstand,scaling(k),rotate(k),trans(:,k),sortIdx(:,k)] = ...
      o.standardizationStep(Xold(:,k),Nnet);
  % Prepare input for advection network
  Xinput(:,1,1,k) = o.prepareInputForNet(Xstand,'advection');
end
% Take a step due to advaction
Xmid = o.translateVinfwNN(Xinput,vback + farFieldtracJump,Xold,rotate,sortIdx);

% Total background velocity
vbackTot = vback + farFieldtracJump;

if ~o.variableKbDt
rotate = zeros(nv,1); trans = zeros(2,nv); sortIdx = zeros(Nnet,nv);
scaling = zeros(nv,1); Xinput = zeros(o.nCompRelax,1,1,nv);
for k = 1 : nv
  % 2) RELAXATION w/ NETWORK
  % Standardize vesicle Xmid
  [Xstand,scaling(k),rotate(k),trans(:,k),sortIdx(:,k)] = ...
      o.standardizationStep(Xmid(:,k),Nnet);
  % Prepare input for relaxation network
  Xinput(:,1,1,k) = o.prepareInputForNet(Xstand,'relaxation');
end
% Take a step
Xnew = o.relaxWNN(Xinput,scaling,rotate,trans,sortIdx,N);
else
Xnew = o.relaxWNNvariableKbDt(Xmid,N,Nnet);    
end

end % DNNsolve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xnew,tenNew,etaNew,RSnew,vbackTot] = newDNNLikeSolve(o,Xold,tenOld,etaOld,RSold,Nnet)
tt = o.tt;
op = tt.op;
if o.confined
  opWall = tt.opWall;
  invWallsMat = tt.bdiagWall;
  walls = o.walls;
  uWalls = walls.u;
  nvbd = walls.nv;
  Nbd = walls.N;
  etaNew = zeros(2*Nbd,nvbd);
  RSnew = zeros(3,nvbd);
else
  vback = o.vinf(Xold);
  etaNew = []; RSnew = [];
end

% semi-implicit time stepping w/ splitting
vesicle = capsules(Xold,[],[],o.kappa,1,1);
vesicle.setUpRate();
nv = vesicle.nv;
N = vesicle.N;
SLPdiag = op.stokesSLmatrixNoCorr(vesicle);
if tt.fmm
  Nup = op.LPuprate * vesicle.N;
  Xup = [interpft(Xold(1:end/2,:),Nup);interpft(Xold(end/2+1:end,:),Nup)];
  vesicleUp = capsules(Xup,[],[],vesicle.kappa,vesicle.viscCont,0); 
  SLPnoCorr = op.stokesSLmatrixNoCorr(vesicleUp);
else
  SLPnoCorr = [];
end
G = op.stokesSLmatrix(vesicle);
[Ben,Ten,Div] = vesicle.computeDerivs;

% Get the near structure, so that we know what to neglect
if o.confined
  [NearV2V,NearV2W] = vesicle.getZone(walls,3);
  [~,NearW2V] = walls.getZone(vesicle,2);
  if nvbd > 1
    NearW2W = o.NearW2W;  
  end % if nvbd > 1
else
  NearV2V = vesicle.getZone([],1);    
end % if o.confined

% Compute bending forces + old tension forces
fBend = vesicle.tracJump(Xold,zeros(N,nv));
fTen = vesicle.tracJump(zeros(2*N,nv),tenOld);
tracJump = fBend+fTen;

% 1) EXPLICIT TENSION AT THE CURRENT STEP
% neglect near interactions

if nv > 1
  % Far-field due to traction jump
  kernelDirect = @op.exactStokesSL;
  if tt.fmm
    kernel = @op.exactStokesSLfmm;
  else
    kernel = @op.exactStokesSL;
  end
  
  
  tt.Galpert = op.stokesSLmatrix(vesicle);  
  SLP = @(X) op.exactStokesSLdiag(vesicle,tt.Galpert,X);
  farFieldtracJump = op.nearSingInt(vesicle,tracJump,SLP,SLPnoCorr,...
      NearV2V,kernel,kernelDirect,vesicle,true,false);
  % if repulsion is on, compute velocity due to that
  if o.repulsion
    if o.confined
      repulsion = vesicle.repulsionSchemeSimple(Xold,o.repStrength,o.repLenScale,...
          walls,[],[]);
    else
      repulsion = vesicle.repulsionSchemeSimple(Xold,o.repStrength,o.repLenScale,...
          [],[],[]);
    end
    Frepulsion = op.exactStokesSLdiag(vesicle,tt.Galpert,repulsion) + ...
        op.nearSingInt(vesicle,repulsion,SLP,SLPnoCorr,...
        NearV2V,kernel,kernelDirect,vesicle,true,false);
    farFieldtracJump = farFieldtracJump + Frepulsion; 
  end % if o.repulsion

else
  farFieldtracJump = zeros(2*N,1);
end

% Far-field due to walls or background flow
if o.confined
  kernelDirect = @opWall.exactStokesDL;
  if tt.fmmDLP
    kernel = @opWall.exactStokesDLnewfmm;
  else
    kernel = @opWall.exactStokesDL;
  end
  DLP = @(X) -1/2*X + opWall.exactStokesDLdiag(walls,tt.wallDLP,X); 
  vback = opWall.nearSingInt(walls,etaOld,DLP,[],NearW2V,kernel,kernelDirect,...
      vesicle,false,false);
  for k = 2:nvbd
    stokeslet = RSold(1:2,k); rotlet = RSold(3,k);  
    vback = vback + tt.RSlets(vesicle.X,walls.center(:,k),stokeslet,rotlet);
  end
end


% SOLVE FOR TENSION
tenNew = zeros(N,nv);
for k = 1 : nv
  LHS = (Div(:,:,k)*G(:,:,k)*Ten(:,:,k));
  selfBend = G(:,:,k)*fBend(:,k);
  RHS = -Div(:,:,k)*(vback(:,k)+farFieldtracJump(:,k)+selfBend);
  tenNew(:,k) = LHS\RHS;
end
% compute force due to tension
fTen = vesicle.tracJump(zeros(2*N,nv),tenNew); tracJump = fBend+fTen;

% 2) SOLVE FOR DENSITY and RS ON WALLS
if o.confined
  % vesicle2wall interactions
  kernelDirect = @op.exactStokesSL;
  if tt.fmm
    kernel = @op.exactStokesSLfmm;
  else
    kernel = @op.exactStokesSL;
  end
  
  SLP = @(X) op.exactStokesSLdiag(vesicle,tt.Galpert,X);  
  ves2wallInt = op.nearSingInt(vesicle,tracJump,SLP,[],NearV2W,...
      kernel,kernelDirect,walls,false,false);
  if o.repulsion
    FrepWall = op.nearSingInt(vesicle,repulsion,SLP,[],NearV2W,kernel,...
        kernelDirect,walls,false,false);
    ves2wallInt = ves2wallInt + FrepWall;
  end % if o.repulsion

  RHS = [uWalls(:)-ves2wallInt(:);zeros(3*(nvbd-1),1)];
  etaRS = invWallsMat*RHS;
  for iw = 1 : nvbd
    etaNew(:,iw) = etaRS((iw-1)*2*Nbd+1:iw*2*Nbd);  
    if iw <= nvbd-1
      RSnew(:,iw+1) = etaRS(2*Nbd*nvbd+(iw-1)*3+1:2*Nbd*nvbd+iw*3);
    end
  end 
  % Update vback due to new density,rotlets and stokeslet
  kernelDirect = @opWall.exactStokesDL;
  if tt.fmmDLP
    kernel = @opWall.exactStokesDLnewfmm;
  else
    kernel = @opWall.exactStokesDL;
  end
  
  
  DLP = @(X) -1/2*X + opWall.exactStokesDLdiag(walls,tt.wallDLP,X);   
  vback = opWall.nearSingInt(walls,etaNew,DLP,[],NearW2V,kernel,...
      kernelDirect,vesicle,false,false);
  
  for k = 2:nvbd
    stokeslet = RSnew(1:2,k); rotlet = RSnew(3,k);  
    vback = vback + tt.RSlets(vesicle.X,walls.center(:,k),stokeslet,rotlet);
  end
end % if iconfined    

% 3) SOLVE FOR NEXT TIME STEP USING OPERATOR SPLITTING
% EXPLICIT-TREATMENT OF INTER-VESICLE BENDING AND TENSION
% TENSION AND X are coupled (semi-implicit).
Xnew = zeros(size(Xold));
if nv > 1
  % compute traction jump due to explicit tension and old bending
  kernelDirect = @op.exactStokesSL;
  if tt.fmm
    kernel = @op.exactStokesSLfmm;
  else
    kernel = @op.exactStokesSL;
  end
  SLP = @(X) op.exactStokesSLdiag(vesicle,tt.Galpert,X);    
  farFieldtracJump = op.nearSingInt(vesicle,tracJump,SLP,SLPnoCorr,NearV2V,kernel,...
      kernelDirect,vesicle,true,false);
  % if repulsion is on, add the velocity (alread computed) due to that
  if o.repulsion
    farFieldtracJump = farFieldtracJump + Frepulsion; 
  end % if o.repulsion
else
  farFieldtracJump = zeros(2*N,1);  
end

% Total background velocity
vbackTot = vback + farFieldtracJump;
% Take a time step
for k = 1 : nv
  M = G(:,:,k)*Ten(:,:,k)*((Div(:,:,k)*G(:,:,k)*Ten(:,:,k))\eye(vesicle.N))*Div(:,:,k);
  rhs1 = Xold(:,k);
  rhs2 = o.dt*(eye(2*vesicle.N)-M)*vbackTot(:,k);
  LHS = (eye(2*vesicle.N)-vesicle.kappa*o.dt*(-G(:,:,k)*Ben(:,:,k)+M*G(:,:,k)*Ben(:,:,k)));
  LHSinv = LHS\eye(2*vesicle.N);
%   Xnew(:,k) = LHS\(rhs1+rhs2);
  Xnew(:,k) = LHSinv * rhs1 + LHSinv*rhs2;
end

end % newDNNLikeSolve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xnew,tenNew,etaNew,RSnew,vbackTot] = DNNsolve2(o,Xold,tenOld,etaOld,RSold,Nnet)
tt = o.tt;
op = tt.op;
if o.confined
  opWall = tt.opWall;
  invWallsMat = tt.bdiagWall;
  walls = o.walls;
  uWalls = walls.u;
  nvbd = walls.nv;
  Nbd = walls.N;
  etaNew = zeros(2*Nbd,nvbd);
  RSnew = zeros(3,nvbd);
else
  vback = o.vinf(Xold);
  etaNew = []; RSnew = [];
end

% FIRST SOLVE FOR Xn+1 and tension
% Then update density


% semi-implicit time stepping w/ splitting
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

% Get the near structure, so that we know what to neglect
if o.confined
  [NearV2V,NearV2W] = vesicle.getZone(walls,3);
  [~,NearW2V] = walls.getZone(vesicle,2);
  if nvbd > 1
    NearW2W = o.NearW2W;  
  end % if nvbd > 1
else
  NearV2V = vesicle.getZone([],1);    
end % if o.confined

% Compute bending forces + old tension forces
tracJump = vesicle.tracJump(Xold,tenOld);
totForce = tracJump;
% 1) EXPLICIT TENSION AT THE CURRENT STEP

if nv > 1
  % Far-field due to traction jump
  kernelDirect = @op.exactStokesSL;
  if tt.fmm
    kernel = @op.exactStokesSLfmm;
  else
    kernel = @op.exactStokesSL;
  end
   
  tt.Galpert = op.stokesSLmatrix(vesicle);  
  SLP = @(X) op.exactStokesSLdiag(vesicle,tt.Galpert,X);
  
  % if repulsion is on, compute velocity due to that
  if o.repulsion
    if o.confined
      repulsion = vesicle.repulsionSchemeSimple(Xold,o.repStrength,o.repLenScale,...
          walls,[],[]);
    else
      repulsion = vesicle.repulsionSchemeSimple(Xold,o.repStrength,o.repLenScale,...
          [],[],[]);
    end
    totForce = totForce + repulsion;
    ves2wallIntRep = op.nearSingInt(vesicle,repulsion,SLP,[],NearV2W,...
        kernel,kernelDirect,walls,false,false);
  end
  farFieldtracJump = op.nearSingInt(vesicle,totForce,SLP,SLPnoCorr,...
      NearV2V,kernel,kernelDirect,vesicle,true,false);
  
  if o.repulsion
    farFieldtracJump = farFieldtracJump + op.exactStokesSLdiag(vesicle,tt.Galpert,repulsion);
  end
else
  farFieldtracJump = zeros(2*N,1);
end

% Far-field due to walls or background flow
if o.confined
  kernelDirect = @opWall.exactStokesDL;
  if tt.fmmDLP
    kernel = @opWall.exactStokesDLnewfmm;
  else
    kernel = @opWall.exactStokesDL;
  end
  DLP = @(X) -1/2*X + opWall.exactStokesDLdiag(walls,tt.wallDLP,X); 
  vback = opWall.nearSingInt(walls,etaOld,DLP,[],NearW2V,kernel,kernelDirect,...
      vesicle,false,false);
  for k = 2:nvbd
    stokeslet = RSold(1:2,k); rotlet = RSold(3,k);  
    vback = vback + tt.RSlets(vesicle.X,walls.center(:,k),stokeslet,rotlet);
  end
end

% farFieldtracJump + vback is the total velocity on a vesicle
% Total background velocity
vbackTot = vback + farFieldtracJump;

% 1) SOLVE FOR NEXT TIME STEP USING OPERATOR SPLITTING
% EXPLICIT-TREATMENT OF INTER-VESICLE BENDING AND TENSION
% TENSION AND X are coupled (semi-implicit).
rotate = zeros(nv,1); trans = zeros(2,nv); sortIdx = zeros(Nnet,nv);
scaling = zeros(nv,1); Xinput = zeros(o.nComp,1,1,nv);
XstandTenExp = []; XinputFourier = zeros(o.nComp,1,1,nv);
for k = 1 : nv
  % 1) TRANSLATION W/ NETWORK
  % Standardize vesicle
  [Xstand,scaling(k),rotate(k),trans(:,k),sortIdx(:,k)] = ...
      o.standardizationStep(Xold(:,k),Nnet);
  
  XinputFourier(:,:,:,k) = o.prepareInputForNet(Xstand,'tensionOnFourier');
  
  % Prepare input for advection network
  Xinput(:,1,1,k) = o.prepareInputForNet(Xstand,'advection');
end
% Take a step due to advaction
Xmid = o.translateVinfwNN(Xinput,vbackTot,Xold,rotate,sortIdx);
scalingTenExp = scaling;
rotateTenExp = rotate;
sortIdxTenExp = sortIdx;

% 2) RELAXATION w/ NETWORK
Xnew = o.relaxWNNvariableKbDt(Xmid,N,Nnet);    



% 2) SOLVE FOR TENSION IN THE NEXT TIME STEP
scaling = zeros(nv,1); rotate = zeros(nv,1); trans = zeros(2,nv);
sortIdx = zeros(Nnet,nv);

XinputTenBen = zeros(o.nComp,1,1,nv);
for k = 1 : nv
  % Standardize vesicle  
  [Xstand,scaling(k),rotate(k),trans(:,k),sortIdx(:,k)] = ...
      o.standardizationStep(Xnew(:,k),Nnet);
  % Prepare inputs for Fourier network and self-bending network
%   XinputFourier(:,:,:,k) = o.prepareInputForNet(Xstand,'tensionOnFourier');
  
  XinputTenBen(:,:,:,k) = o.prepareInputForNet(Xstand,'tensionSelfBend');
end
% Approximate inv(Div*G*Ten)*Div*vExt
% vBackSolve = o.invTenMatOnVback(XinputFourier,vbackTot,scaling,rotate,sortIdx);
vBackSolve = o.invTenMatOnVback(XinputFourier,vbackTot,scalingTenExp,rotateTenExp,sortIdxTenExp);
% Approximate inv(Div*G*Ten)*G*(-Ben)*x
selfBendSolve = o.invTenMatOnSelfBend(XinputTenBen,scaling,sortIdx,N);
% update tension
tenNew = -(vBackSolve+selfBendSolve);


% Update the velocities on the walls and solve for the walls

vesicle = capsules(Xnew,[],[],o.kappa,1,1);
vesicle.setUpRate();
% 
% % Get the near structure, so that we know what to neglect
if o.confined
  [~,NearV2W] = vesicle.getZone(walls,3);
end % if o.confined

% Update the traction jump
tracJump = vesicle.tracJump(Xnew,tenNew); 

% Update tt.Galpert, SLP
tt.Galpert = op.stokesSLmatrix(vesicle);  

if o.confined
  % 3) SOLVE FOR DENSITY and RS ON WALLS
  % vesicle2wall interactions
  kernelDirect = @op.exactStokesSL;
  if tt.fmm
    kernel = @op.exactStokesSLfmm;
  else
    kernel = @op.exactStokesSL;
  end

  SLP = @(X) op.exactStokesSLdiag(vesicle,tt.Galpert,X); 
  ves2wallInt = op.nearSingInt(vesicle,tracJump,SLP,[],NearV2W,...
        kernel,kernelDirect,walls,false,false);
  if o.repulsion
  RHS = [uWalls(:)-ves2wallInt(:)-ves2wallIntRep(:);zeros(3*(nvbd-1),1)];
  else
  RHS = [uWalls(:)-ves2wallInt(:);zeros(3*(nvbd-1),1)];
  end
  etaRS = invWallsMat*RHS;
  for iw = 1 : nvbd
    etaNew(:,iw) = etaRS((iw-1)*2*Nbd+1:iw*2*Nbd);  
    if iw <= nvbd-1
      RSnew(:,iw+1) = etaRS(2*Nbd*nvbd+(iw-1)*3+1:2*Nbd*nvbd+iw*3);
    end
  end 
end % if iconfined    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vBackSolve = invTenMatOnVback(o,Xinput,vext,scaling,rotate,sortIdx)
% Approximate inv(Div*G*Ten)*Div*vExt 
    
% number of vesicles
nv = numel(vext(1,:));
% number of points of exact solve
N = numel(vext(:,1))/2;    
% number of points for network
Nnet = numel(sortIdx(:,1));
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
Nnet = numel(sortIdx(:,1));
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
if o.nCompRelax > 32
  % predict PCA coefficients for Xnew
  Ypred(:,33:48) = predict(nets{3},Xinput(33:48,:,:,:));
  % we use normalized output so take that back
  Ypred(:,33:48) = (Ypred(:,33:48)-offset(3))*sdevChan1(3)/scale(3)+muChan1(3); 
end
if o.nCompRelax > 48
  % predict PCA coefficients for Xnew
  Ypred(:,49:64) = predict(nets{4},Xinput(49:64,:,:,:));
  % we use normalized output so take that back
  Ypred(:,49:64) = (Ypred(:,49:64)-offset(4))*sdevChan1(4)/scale(4)+muChan1(4); 
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
function Xnew = relaxWNNvariableKbDt(o,Xmid,N,Nnet)
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
      o.standardizationStep(Xmid(:,k),Nnet);
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
function [X,scaling,rotate,trans,sortIdx] = standardizationStep(o,Xin,Nnet)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xnew,tenNew,etaNew,RSnew] = DNNlikeExactSolve(o,Xold,tenOld,etaOld,RSold)
 
tt = o.tt;
op = tt.op;
if o.confined
  opWall = tt.opWall;
  invWallsMat = tt.bdiagWall;
  walls = o.walls;
  uWalls = walls.u;
  nvbd = walls.nv;
  Nbd = walls.N;
  etaNew = zeros(2*Nbd,nvbd);
  RSnew = zeros(3,nvbd);
else
  vback = o.vinf(Xold);
  etaNew = []; RSnew = [];
end

% Semi-implicit time stepping w/splitting
vesicle = capsules(Xold,[],[],o.kappa,1,1);
vesicle.setUpRate();
nv = vesicle.nv;
N = vesicle.N;
% Build the operators
G = op.stokesSLmatrix(vesicle);
[~,Ten,Div] = vesicle.computeDerivs;
if tt.fmm
  Nup = op.LPuprate * vesicle.N;
  Xup = [interpft(Xold(1:end/2,:),Nup);interpft(Xold(end/2+1:end,:),Nup)];
  vesicleUp = capsules(Xup,[],[],vesicle.kappa,vesicle.viscCont,0); 
  SLPnoCorr = op.stokesSLmatrixNoCorr(vesicleUp);
else
  SLPnoCorr = [];
end

% Compute bending forces + old tension forces
fBend = vesicle.tracJump(Xold,zeros(N,nv));
fTen = vesicle.tracJump(zeros(2*N,nv),tenOld);
tracJump = fBend+fTen;

% 1) EXPLICIT TENSION AT THE CURRENT STEP
% neglect near interactions

% Get the near structure, so that we know what to neglect
if o.confined
  [NearV2V,NearV2W] = vesicle.getZone(walls,3);
  [~,NearW2V] = walls.getZone(vesicle,2);
  if nvbd > 1
    NearW2W = o.NearW2W;  
  end % if nvbd > 1
else
  NearV2V = vesicle.getZone([],1);    
end % if o.confined

if nv > 1
% Far-field due to traction jump
kernelDirect = @op.exactStokesSL;
if tt.fmm
  kernel = @op.exactStokesSL;
else
  kernel = @op.exactStokesSL;
end

if o.useTrueNear
  farFieldtracJump = op.nearSingInt(vesicle,tracJump,SLPnoCorr,NearV2V,...
      kernel,kernelDirect,vesicle,true,false);
else
  [~,nearField,farFieldtracJump] = op.divideNearFarSLP(vesicle,tracJump,SLPnoCorr,...
      NearV2V,kernel,kernelDirect,vesicle,true);      
  if o.useNear
    farFieldtracJump = farFieldtracJump + nearField;
  end
end

else
farFieldtracJump = zeros(2*N,1);
end

% Far-field due to walls or background flow
if o.confined
  kernelDirect = @opWall.exactStokesDL;
  if tt.fmmDLP
    kernel = @opWall.exactStokesDL;
  else
    kernel = @opWall.exactStokesDL;
  end
  if o.useTrueNear
    vback = opWall.nearSingInt(walls,etaOld,[],NearW2V,kernel,...
        kernelDirect,vesicle,false,false);
  else
    [~,nearField,vback] = opWall.divideNearFarSLP(walls,etaOld,[],NearW2V,...
        kernel,kernelDirect,vesicle,false);
    if o.useNear
      vback = vback + nearField;
    end
  end
  
  for k = 2:nvbd
    stokeslet = RSold(1:2,k); rotlet = RSold(3,k);  
    vback = vback + tt.RSlets(vesicle.X,walls.center(:,k),stokeslet,rotlet);
  end
end

tenNew = zeros(N,nv);
for k = 1 : nv
  LHS = (Div(:,:,k)*G(:,:,k)*Ten(:,:,k));
  % DNN approximates inv(Div*G*Ten)*G*(-Ben)*x and inv(Div*G*Ten)*Div*Fourier
  
  % SLP due to bending force on itself 
  selfBend = G(:,:,k)*fBend(:,k);
  
  RHS = -Div(:,:,k)*(vback(:,k)+farFieldtracJump(:,k)+selfBend);
  %[inv(LHS)*(Div(:,:,k)*(vback(:,k)+farFieldtracJump(:,k))) inv(LHS)*(Div(:,:,k)*selfBend)]
  %pause
  tenNew(:,k) = LHS\RHS;
end
% compute force due to tension
fTen = vesicle.tracJump(zeros(2*N,nv),tenNew); tracJump = fBend+fTen;

if o.confined
  % 2) SOLVE FOR DENSITY and RS ON WALLS
  
  % vesicle2wall interactions
  kernelDirect = @op.exactStokesSL;
  if tt.fmm
    kernel = @op.exactStokesSLfmm;
  else
    kernel = @op.exactStokesSL;
  end
  [~,nearField,ves2wallInt] = op.divideNearFarSLP(vesicle,tracJump,[],NearV2W,...
      kernel,kernelDirect,walls,0);
  if o.useNear
    ves2wallInt = ves2wallInt + nearField;
  end
  RHS = [uWalls(:)-ves2wallInt(:);zeros(3*(nvbd-1),1)];
  etaRS = invWallsMat*RHS;
  for iw = 1 : nvbd
    etaNew(:,iw) = etaRS((iw-1)*2*Nbd+1:iw*2*Nbd);  
    if iw <= nvbd-1
      RSnew(:,iw+1) = etaRS(2*Nbd*nvbd+(iw-1)*3+1:2*Nbd*nvbd+iw*3);
    end
  end
  % Update vback due to walls
  kernelDirect = @opWall.exactStokesDL;
  if tt.fmmDLP
    kernel = @opWall.exactStokesDL;
  else
    kernel = @opWall.exactStokesDL;
  end
  [~,nearField,vback] = opWall.divideNearFarSLP(walls,etaNew,[],NearW2V,...
      kernel,kernelDirect,vesicle,false);
  if o.useNear
    vback = vback + nearField;
  end
  for k = 2:nvbd
    stokeslet = RSnew(1:2,k); rotlet = RSnew(3,k);  
    vback = vback + tt.RSlets(vesicle.X,walls.center(:,k),stokeslet,rotlet);
  end
end % if iconfined

% 3) SOLVE FOR NEXT TIME STEP USING OPERATOR SPLITTING
% EXPLICIT-TREATMENT OF INTER-VESICLE BENDING AND TENSION
% TENSION AND X are coupled (semi-implicit).
Xnew = zeros(size(Xold));
if nv > 1
% compute traction jump due to explicit tension and old bending
kernelDirect = @op.exactStokesSL;
if tt.fmm
  kernel = @op.exactStokesSL;
else
  kernel = @op.exactStokesSL;
end
[~,nearField,farFieldtracJump] = op.divideNearFarSLP(...
        vesicle,tracJump,[],NearV2V,kernel,kernelDirect,vesicle,1);
if o.useNear 
  farFieldtracJump = farFieldtracJump + nearField;
end
else
farFieldtracJump = zeros(2*N,1);
end

% Walls2Ves is already computed if confined (vback)    
for k = 1 : nv
  % First, translate with vext
  vext = vback(:,k) + farFieldtracJump(:,k);
  %Xmid = o.translateVinf(vext,o.dt,Xold(:,k),op);
  
  % Second, solve bending problem
  vesicle = capsules(Xold(:,k),[],[],1,1,1); vesicle.setUpRate();
  Xnew(:,k) = o.relaxExactSolve(vesicle,vext,o.dt,Xold(:,k),op);
end


end % DNNlikeExactSolve
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
end
    
end % setBgFlow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tt = buildTstep(o,X,pramsIn)
N = pramsIn.N; nv = pramsIn.nv;

prams.N = N;
prams.nv = nv;
prams.Nbd = pramsIn.Nbd;
prams.nvbd = pramsIn.nvbd;
prams.kappa = pramsIn.kappa;
options.diffDiscWalls = 0;
prams.NbdInt = 0; prams.NbdExt = 0; prams.nvbdInt = 0; prams.nvbdExt = 0;
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

options.correctShape = 1;
options.adhesion = 0;
options.repulsion = 0;
if strcmp(pramsIn.bgFlow,'couette')
options.confined = true;
options.farField = 'couette';
options.farFieldSpeed = pramsIn.speed;
else
options.confined = false;
end

[options,prams] = initVes2D(options,prams);
    
om = monitor(X,options,prams);
tt = tstep(options,prams,om);
if ~options.confined % if not confined, then assign vinf
tt.farField = o.setBgFlow(pramsIn.bgFlow,pramsIn.speed); 
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % methods

end % dnnTools
