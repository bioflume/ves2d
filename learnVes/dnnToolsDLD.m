classdef dnnToolsDLD

properties
KbDts
interpOrder
repLenScale
variableKbDt
ttExact
Nexact
N
confined
tt
dt
vinf
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
oc
useNear
useTrueNear
repulsion
repStrength
minDist
kappa
opRigid
DLPEll
invTenMatEll
selfBendMatEll
Mell
relaxMatEll
Nrigid

end

methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = dnnToolsDLD(X,XwallsInt,XwallsExt,prams)

o.confined = true;
o.kappa = prams.kappa;
o.dt = prams.dt;    
o.interpOrder = prams.interpOrder;
o.Nexact = prams.Nexact;
o.N = prams.N;
o.tt = o.buildTstep(X,prams);  
prams.N = o.Nexact;
o.ttExact = o.buildTstep(X,prams);
prams.N = o.N;
o = o.initWalls(XwallsInt,XwallsExt,prams);

end % dnnToolsDLD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = initWalls(o,XwallsInt,XwallsExt,prams)
tt = o.tt;
% poten classes for walls
potWallInt = tt.opWallInt;
potWallExt = tt.opWallExt;
% velocity on solid walls coming from no slip boundary condition
[uwallsExt,uwallsInt] = tt.farField(XwallsExt,XwallsInt);


% build walls
wallsInt = capsules(XwallsInt,[],uwallsInt,...
      zeros(prams.nvbdInt,1),zeros(prams.nvbdInt,1),1);
wallsExt = capsules(XwallsExt,[],uwallsExt,0,0,1);
wallsInt.setUpRate(potWallInt);
wallsExt.setUpRate(potWallExt);
o.wallsInt = wallsInt; o.wallsExt = wallsExt;

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

ttExact = o.ttExact;
if isempty(ttExact.wallDLPint)
  ttExact.wallDLPint = potWallInt.stokesDLmatrix(wallsInt);
end
if isempty(ttExact.wallDLPext)
  ttExact.wallDLPext = potWallExt.stokesDLmatrix(wallsExt);
end

if isempty(ttExact.wallDLPintNoCorr)
  ttExact.wallDLPintNoCorr = potWallInt.stokesDLmatrixNoCorr(wallsInt);
end
if isempty(ttExact.wallDLPextNoCorr)
  ttExact.wallDLPextNoCorr = potWallExt.stokesDLmatrixNoCorr(wallsExt);    
end


% N0 to remove rank-1 deficiency
if isempty(ttExact.wallN0)
  ttExact.wallN0 = potWallExt.stokesN0matrix(wallsExt);
end


% inverse of the walls matrix (-1/2*density + DLP + N0 + R-S)
tt.matFreeWalls = true; % so, create walls matrix and save it
tt.bdiagWall = tt.wallsPrecondDiffDisc(wallsInt,wallsExt);
ttExact.matFreeWalls = tt.matFreeWalls;
ttExact.bdiagWall = tt.bdiagWall;
o.ttExact = ttExact;
o.tt = tt; 
% tt includes tt.wallDLPandRSmat which is walls matrix saved 

end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xnew,tenNew,etaIntNew,etaExtNew,RSnew,vbackTot] = DNNsolveMixed(o,Xold,...
        tenOld,etaIntOld,etaExtOld,RSold,Nnet)
    
tt = o.tt;
op = tt.op;

opWallInt = tt.opWallInt;
opWallExt = tt.opWallExt;  
wallsInt = o.wallsInt;
wallsExt = o.wallsExt;
nvbdInt = wallsInt.nv;
nvbdExt = 1;
nvbd = nvbdExt + nvbdInt;
NbdInt = wallsInt.N;
NbdExt = wallsExt.N;
etaIntNew = zeros(2*NbdInt,nvbdInt);
etaExtNew = zeros(2*NbdExt,nvbdExt);
RSnew = zeros(3,nvbd);

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


% Compute bending forces + old tension forces
fBend = vesicle.tracJump(Xold,zeros(N,nv));
fTen = vesicle.tracJump(zeros(2*N,nv),tenOld);
tracJump = fBend+fTen;

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
  farFieldtracJump = op.nearSingInt(vesicle,tracJump,SLP,SLPnoCorr,...
      NearV2V,kernel,kernelDirect,vesicle,true,false);
else
  farFieldtracJump = zeros(2*N,1);   
end

if o.repulsion
  repulsion = vesicle.repulsionScheme(Xold,o.repStrength,o.minDist,[],wallsInt,wallsExt);
  Frepulsion = op.exactStokesSLdiag(vesicle,tt.Galpert,repulsion) + ...
      op.nearSingInt(vesicle,repulsion,SLP,SLPnoCorr,NearV2V,kernel,kernelDirect,vesicle,true,true);
  
  FREPwallInt = op.nearSingInt(vesicle,repulsion,SLP,[],NearV2Wint,kernel,kernelDirect,wallsInt,false,false);
  FREPwallExt = op.nearSingInt(vesicle,repulsion,SLP,[],NearV2Wext,kernel,kernelDirect,wallsExt,false,false);
else
  Frepulsion = zeros(2*N,nv);  
  FREPwallExt = zeros(2*NbdExt,nvbdExt);
  FREPwallInt = zeros(2*NbdInt,nvbdInt);
end % if o.repulsion

% Far-field due to walls or background flow
if o.confined
  kernelDirect = @opWallInt.exactStokesDL;
  kernelDirect2 = @opWallExt.exactStokesDL;
  if tt.fmmDLP
    kernel = @opWallInt.exactStokesDLnewfmm;
    kernel2 = @opWallExt.exactStokesDLnewfmm;
  else
    kernel = @opWallInt.exactStokesDL;
    kernel2 = @opWallExt.exactStokesDL;
  end
  DLP = @(X) -1/2*opWallInt.exactStokesDLdiag(wallsInt,tt.wallDLPint,X); 
  vback = opWallInt.nearSingInt(wallsInt,etaIntOld,DLP,[],NearWint2V,...
      kernel,kernelDirect,vesicle,false,false);
  
  DLP = @(X) -1/2*opWallExt.exactStokesDLdiag(wallsExt,tt.wallDLPext,X);
  vbackExt = opWallExt.nearSingInt(wallsExt,etaExtOld,DLP,[],NearWext2V,...
      kernel2,kernelDirect2,vesicle,false,false);
  vback = vback + vbackExt;
  if o.useNear
    vback = vback + nearField;
  end

  for k = 2:nvbd
    stokeslet = RSold(1:2,k); rotlet = RSold(3,k);  
    vback = vback + tt.RSlets(vesicle.X,wallsInt.center(:,k-1),stokeslet,rotlet);
  end
end

scaling = zeros(nv-1,1); rotate = zeros(nv-1,1); trans = zeros(2,nv-1);
sortIdx = zeros(Nnet,nv-1);
XinputFourier = zeros(o.nComp,1,1,nv-1);
XinputTenBen = zeros(o.nComp,1,1,nv-1);
for k = 2 : nv
  % Standardize vesicle  
  [Xstand,scaling(k-1),rotate(k-1),trans(:,k-1),sortIdx(:,k-1)] = ...
      o.standardizationStep(Xold(:,k),Nnet);
  
  % Prepare inputs for Fourier network and self-bending network
  XinputFourier(:,:,:,k-1) = o.prepareInputForNet(Xstand,'tensionOnFourier');
  XinputTenBen(:,:,:,k-1) = o.prepareInputForNet(Xstand,'tensionSelfBend');

end
% Approximate inv(Div*G*Ten)*Div*vExt
vBackSolve = o.invTenMatOnVback(XinputFourier,vback(:,2:end)+farFieldtracJump(:,2:end),...
    scaling,rotate,sortIdx);
% Approximate inv(Div*G*Ten)*G*(-Ben)*x
selfBendSolve = o.invTenMatOnSelfBend(XinputTenBen,scaling,sortIdx,N);
% update tension
tenNew = -(vBackSolve+selfBendSolve);

% COMPUTE TENSION FOR RA = 0.9
[XEllStand,scalingEll,rotateEll,~,sortIdxEll] = ...
    o.standardizationStep(Xold(:,1),Nnet);      
vext = vback(:,1) + farFieldtracJump(:,1);
vext = [interpft(vext(1:end/2),Nnet);interpft(vext(end/2+1:end),Nnet)];
vextStand = o.standardize(vext,[0;0],rotateEll,1,sortIdxEll);
vBackEllStand = o.invTenMatEll*vextStand;
vBackEll = zeros(size(vBackEllStand));
vBackEll(sortIdxEll) = vBackEllStand;

selfBendStand = o.selfBendMatEll*XEllStand;
selfBendEll = zeros(size(selfBendStand));
selfBendEll(sortIdxEll) = selfBendStand/scalingEll^2;
tenEll = -(vBackEll+selfBendEll);
tenEll = interpft(tenEll,N);
tenNew = [tenEll tenNew];


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
  
  SLP = @(X) op.exactStokesSLdiag(vesicle,tt.Galpert,X);  
  ves2wallInt = op.nearSingInt(vesicle,tracJump,SLP,[],NearV2Wint,...
      kernel,kernelDirect,wallsInt,false,false);
  ves2wallExt = op.nearSingInt(vesicle,tracJump,SLP,[],NearV2Wext,...
      kernel,kernelDirect,wallsExt,false,false);
  
  RHSint = wallsInt.u(:)-ves2wallInt(:)-FREPwallInt(:);
  RHSext = wallsExt.u(:)-ves2wallExt(:)-FREPwallExt(:);
  RHS = [RHSext;RHSint;zeros(3*(nvbd-1),1)];
  etaRS = tt.bdiagWall*RHS;
  etaExtNew = etaRS(1:2*NbdExt);
  etaRS = etaRS(2*NbdExt+1:end);
  for iw = 1 : nvbdInt
    etaIntNew(:,iw) = etaRS((iw-1)*2*NbdInt+1:iw*2*NbdInt);  
    if iw <= nvbdInt
      RSnew(:,iw+1) = etaRS(2*NbdInt*nvbdInt+(iw-1)*3+1:2*NbdInt*nvbdInt+iw*3);
    end
  end 
  % Update vback due to new density,rotlets and stokeslet
  kernelDirect = @opWallInt.exactStokesDL;
  kernelDirect2 = @opWallExt.exactStokesDL;
  if tt.fmmDLP
    kernel = @opWallInt.exactStokesDLnewfmm;
    kernel2 = @opWallExt.exactStokesDLnewfmm;
  else
    kernel = @opWallInt.exactStokesDL;
    kernel2 = @opWallExt.exactStokesDL;
  end
  
  DLP = @(X) -1/2*X + opWallInt.exactStokesDLdiag(wallsInt,tt.wallDLPint,X);
  vback = opWallInt.nearSingInt(wallsInt,etaIntNew,DLP,[],NearWint2V,...
      kernel,kernelDirect,vesicle,false,false);
  
  DLP = @(X) -1/2*X + opWallExt.exactStokesDLdiag(wallsExt,tt.wallDLPext,X);
  vbackExt = opWallExt.nearSingInt(wallsExt,etaExtNew,DLP,[],NearWext2V,...
      kernel2,kernelDirect2,vesicle,false);
  vback = vback + vbackExt;
  
  for k = 2:nvbd
    stokeslet = RSnew(1:2,k); rotlet = RSnew(3,k);  
    vback = vback + tt.RSlets(vesicle.X,wallsInt.center(:,k-1),stokeslet,rotlet);
  end
end % if iconfined    

% 3) SOLVE FOR NEXT TIME STEP USING OPERATOR SPLITTING
% EXPLICIT-TREATMENT OF INTER-VESICLE BENDING AND TENSION
% TENSION AND X are coupled (semi-implicit).
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
else
  farFieldtracJump = zeros(2*N,1);  
end
farFieldtracJump = farFieldtracJump + Frepulsion;

% Walls2Ves is already computed if confined (vback)    
rotate = zeros(nv-1,1); trans = zeros(2,nv-1); sortIdx = zeros(Nnet,nv-1);
scaling = zeros(nv-1,1); Xinput = zeros(o.nComp,1,1,nv-1);
for k = 2 : nv
  % 1) TRANSLATION W/ NETWORK
  % Standardize vesicle
  [Xstand,scaling(k-1),rotate(k-1),trans(:,k-1),sortIdx(:,k-1)] = ...
      o.standardizationStep(Xold(:,k),Nnet);
  % Prepare input for advection network
  Xinput(:,1,1,k-1) = o.prepareInputForNet(Xstand,'advection');
end
% Take a step due to advaction
Xmid = o.translateVinfwNN(Xinput,vback(:,2:end) + farFieldtracJump(:,2:end),...
  Xold(:,2:end),rotate,sortIdx);

vbackTot = vback + farFieldtracJump;

% Compute for RA = 0.9
[~,~,rotateEll,~,sortIdxEll] = ...
    o.standardizationStep(Xold(:,1),Nnet);
vext = vback(:,1)+farFieldtracJump(:,1);
vextUp = [interpft(vext(1:end/2),Nnet);interpft(vext(end/2+1:end),Nnet)];
vextStand = o.standardize(vextUp,[0;0],rotateEll,1,sortIdxEll);
MVinfEllStand = o.Mell*vextStand;
MVinfEll = zeros(size(MVinfEllStand));
MVinfEll([sortIdxEll;sortIdxEll+Nnet]) = MVinfEllStand;
MVinfEll = o.rotationOperator(MVinfEll,-rotateEll);
MVinfEll = [interpft(MVinfEll(1:end/2),N);interpft(MVinfEll(end/2+1:end),N)];
XmidEll = Xold(:,1) + o.dt*vext - o.dt*MVinfEll;



rotate = zeros(nv-1,1); trans = zeros(2,nv-1); sortIdx = zeros(Nnet,nv-1);
scaling = zeros(nv-1,1); Xinput = zeros(o.nCompRelax,1,1,nv-1);
for k = 2 : nv
  % 2) RELAXATION w/ NETWORK
  % Standardize vesicle Xmid
  [Xstand,scaling(k-1),rotate(k-1),trans(:,k-1),sortIdx(:,k-1)] = ...
      o.standardizationStep(Xmid(:,k-1),Nnet);
  % Prepare input for relaxation network
  Xinput(:,1,1,k-1) = o.prepareInputForNet(Xstand,'relaxation');
end
% Take a step
Xnew = o.relaxWNN(Xinput,scaling,rotate,trans,sortIdx,N);

% Compute for RA = 0.9
[XstandEll,scalingEll,rotateEll,transEll,sortIdxEll] = ...
    o.standardizationStep(XmidEll,Nnet);
XnewEll = o.relaxMatEll*XstandEll;
XnewEll = o.destandardize(XnewEll,transEll,rotateEll,scalingEll,sortIdxEll);
XnewEll = [interpft(XnewEll(1:end/2),N);interpft(XnewEll(end/2+1:end),N)];
Xnew = [XnewEll Xnew];
    
end % DNNsolveMixed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xnew,tenNew,etaIntNew,etaExtNew,RSnew,vbackTot] = DNNsolveMixed2(o,Xold,...
        tenOld,etaIntOld,etaExtOld,RSold,Nnet)
    
tt = o.tt;
op = tt.op;

opWallInt = tt.opWallInt;
opWallExt = tt.opWallExt;  
wallsInt = o.wallsInt;
wallsExt = o.wallsExt;
nvbdInt = wallsInt.nv;
nvbdExt = 1;
nvbd = nvbdExt + nvbdInt;
NbdInt = wallsInt.N;
NbdExt = wallsExt.N;
etaIntNew = zeros(2*NbdInt,nvbdInt);
etaExtNew = zeros(2*NbdExt,nvbdExt);
RSnew = zeros(3,nvbd);

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


% Compute bending forces + old tension forces
fBend = vesicle.tracJump(Xold,zeros(N,nv));
fTen = vesicle.tracJump(zeros(2*N,nv),tenOld);
tracJump = fBend+fTen;

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
  farFieldtracJump = op.nearSingInt(vesicle,tracJump,SLP,SLPnoCorr,...
      NearV2V,kernel,kernelDirect,vesicle,true,false);
else
  farFieldtracJump = zeros(2*N,1);   
end

if o.repulsion
  repulsion = vesicle.repulsionSchemeSimple(Xold,o.repStrength,o.repLenScale,[],wallsInt,wallsExt);
  Frepulsion = op.exactStokesSLdiag(vesicle,tt.Galpert,repulsion) + ...
      op.nearSingInt(vesicle,repulsion,SLP,SLPnoCorr,NearV2V,kernel,kernelDirect,vesicle,true,true);
  
  FREPwallInt = op.nearSingInt(vesicle,repulsion,SLP,[],NearV2Wint,kernel,kernelDirect,wallsInt,false,false);
  FREPwallExt = op.nearSingInt(vesicle,repulsion,SLP,[],NearV2Wext,kernel,kernelDirect,wallsExt,false,false);
else
  Frepulsion = zeros(2*N,nv);  
  FREPwallExt = zeros(2*NbdExt,nvbdExt);
  FREPwallInt = zeros(2*NbdInt,nvbdInt);
end % if o.repulsion

% Far-field due to walls or background flow
if o.confined
  kernelDirect = @opWallInt.exactStokesDL;
  kernelDirect2 = @opWallExt.exactStokesDL;
  if tt.fmmDLP
    kernel = @opWallInt.exactStokesDLnewfmm;
    kernel2 = @opWallExt.exactStokesDLnewfmm;
  else
    kernel = @opWallInt.exactStokesDL;
    kernel2 = @opWallExt.exactStokesDL;
  end
  DLP = @(X) -1/2*opWallInt.exactStokesDLdiag(wallsInt,tt.wallDLPint,X); 
  vback = opWallInt.nearSingInt(wallsInt,etaIntOld,DLP,[],NearWint2V,...
      kernel,kernelDirect,vesicle,false,false);
  
  DLP = @(X) -1/2*opWallExt.exactStokesDLdiag(wallsExt,tt.wallDLPext,X);
  vbackExt = opWallExt.nearSingInt(wallsExt,etaExtOld,DLP,[],NearWext2V,...
      kernel2,kernelDirect2,vesicle,false,false);
  vback = vback + vbackExt;
  
  for k = 2:nvbd
    stokeslet = RSold(1:2,k); rotlet = RSold(3,k);  
    vback = vback + tt.RSlets(vesicle.X,wallsInt.center(:,k-1),stokeslet,rotlet);
  end
end

vbackTot = vback + farFieldtracJump;
% Walls2Ves is already computed if confined (vback)    
rotate = zeros(nv-1,1); trans = zeros(2,nv-1); sortIdx = zeros(Nnet,nv-1);
scaling = zeros(nv-1,1); Xinput = zeros(o.nComp,1,1,nv-1);
for k = 2 : nv
  % 1) TRANSLATION W/ NETWORK
  % Standardize vesicle
  [Xstand,scaling(k-1),rotate(k-1),trans(:,k-1),sortIdx(:,k-1)] = ...
      o.standardizationStep(Xold(:,k),Nnet);
  % Prepare input for advection network
  Xinput(:,1,1,k-1) = o.prepareInputForNet(Xstand,'advection');
end
% Take a step due to advaction
Xmid = o.translateVinfwNN(Xinput,vbackTot(:,2:end),...
  Xold(:,2:end),rotate,sortIdx);


% Compute for RA = 0.9
[~,~,rotateEll,~,sortIdxEll] = ...
    o.standardizationStep(Xold(:,1),Nnet);
vext = vbackTot(:,1);
vextUp = [interpft(vext(1:end/2),Nnet);interpft(vext(end/2+1:end),Nnet)];
vextStand = o.standardize(vextUp,[0;0],rotateEll,1,sortIdxEll);
MVinfEllStand = o.Mell*vextStand;
MVinfEll = zeros(size(MVinfEllStand));
MVinfEll([sortIdxEll;sortIdxEll+Nnet]) = MVinfEllStand;
MVinfEll = o.rotationOperator(MVinfEll,-rotateEll);
MVinfEll = [interpft(MVinfEll(1:end/2),N);interpft(MVinfEll(end/2+1:end),N)];
XmidEll = Xold(:,1) + o.dt*vext - o.dt*MVinfEll;

% variable KbDt for Xnew;
Xnew = o.relaxWNNvariableKbDt(Xmid,N,Nnet);


% Compute for RA = 0.9
[XstandEll,scalingEll,rotateEll,transEll,sortIdxEll] = ...
    o.standardizationStep(XmidEll,Nnet);
XnewEll = o.relaxMatEll*XstandEll;
XnewEll = o.destandardize(XnewEll,transEll,rotateEll,scalingEll,sortIdxEll);
XnewEll = [interpft(XnewEll(1:end/2),N);interpft(XnewEll(end/2+1:end),N)];
Xnew = [XnewEll Xnew];


scaling = zeros(nv-1,1); rotate = zeros(nv-1,1); trans = zeros(2,nv-1);
sortIdx = zeros(Nnet,nv-1);
XinputFourier = zeros(o.nComp,1,1,nv-1);
XinputTenBen = zeros(o.nComp,1,1,nv-1);
for k = 2 : nv
  % Standardize vesicle  
  [Xstand,scaling(k-1),rotate(k-1),trans(:,k-1),sortIdx(:,k-1)] = ...
      o.standardizationStep(Xnew(:,k),Nnet);
  
  % Prepare inputs for Fourier network and self-bending network
  XinputFourier(:,:,:,k-1) = o.prepareInputForNet(Xstand,'tensionOnFourier');
  XinputTenBen(:,:,:,k-1) = o.prepareInputForNet(Xstand,'tensionSelfBend');

end
% Approximate inv(Div*G*Ten)*Div*vExt
vBackSolve = o.invTenMatOnVback(XinputFourier,vbackTot(:,2:end),scaling,rotate,sortIdx);
% Approximate inv(Div*G*Ten)*G*(-Ben)*x
selfBendSolve = o.invTenMatOnSelfBend(XinputTenBen,scaling,sortIdx,N);
% update tension
tenNew = -(vBackSolve+selfBendSolve);

% COMPUTE TENSION FOR RA = 0.9
[XEllStand,scalingEll,rotateEll,~,sortIdxEll] = ...
    o.standardizationStep(Xnew(:,1),Nnet);      
vext = vbackTot(:,1);
vext = [interpft(vext(1:end/2),Nnet);interpft(vext(end/2+1:end),Nnet)];
vextStand = o.standardize(vext,[0;0],rotateEll,1,sortIdxEll);
vBackEllStand = o.invTenMatEll*vextStand;
vBackEll = zeros(size(vBackEllStand));
vBackEll(sortIdxEll) = vBackEllStand;

selfBendStand = o.selfBendMatEll*XEllStand;
selfBendEll = zeros(size(selfBendStand));
selfBendEll(sortIdxEll) = selfBendStand/scalingEll^2;
tenEll = -(vBackEll+selfBendEll);
tenEll = interpft(tenEll,N);
tenNew = [tenEll tenNew];


% Update the traction jump
vesicle = capsules(Xnew, [], [], o.kappa,1,1);
vesicle.setUpRate();
tracJump = vesicle.tracJump(Xnew,tenNew);

% Get the near structure, so that we know what to neglect
[~,NearV2Wint] = vesicle.getZone(wallsInt,3);
[~,NearV2Wext] = vesicle.getZone(wallsExt,2);



if o.confined
  % 2) SOLVE FOR DENSITY and RS ON WALLS
  % vesicle2wall interactions
  kernelDirect = @op.exactStokesSL;
  if tt.fmm
    kernel = @op.exactStokesSLfmm;
  else
    kernel = @op.exactStokesSL;
  end
  
  SLP = @(X) op.exactStokesSLdiag(vesicle,tt.Galpert,X);
  
  ves2wallInt = op.nearSingInt(vesicle,tracJump,SLP,[],NearV2Wint,...
      kernel,kernelDirect,wallsInt,false,false);
  ves2wallExt = op.nearSingInt(vesicle,tracJump,SLP,[],NearV2Wext,...
      kernel,kernelDirect,wallsExt,false,false);
  
  RHSint = wallsInt.u(:)-ves2wallInt(:)-FREPwallInt(:);
  RHSext = wallsExt.u(:)-ves2wallExt(:)-FREPwallExt(:);
  RHS = [RHSext;RHSint;zeros(3*(nvbd-1),1)];
  etaRS = tt.bdiagWall*RHS;
  etaExtNew = etaRS(1:2*NbdExt);
  etaRS = etaRS(2*NbdExt+1:end);
  for iw = 1 : nvbdInt
    etaIntNew(:,iw) = etaRS((iw-1)*2*NbdInt+1:iw*2*NbdInt);  
    if iw <= nvbdInt
      RSnew(:,iw+1) = etaRS(2*NbdInt*nvbdInt+(iw-1)*3+1:2*NbdInt*nvbdInt+iw*3);
    end
  end 
  % Update vback due to new density,rotlets and stokeslet
  kernelDirect = @opWallInt.exactStokesDL;
  kernelDirect2 = @opWallExt.exactStokesDL;
  if tt.fmmDLP
    kernel = @opWallInt.exactStokesDLnewfmm;
    kernel2 = @opWallExt.exactStokesDLnewfmm;
  else
    kernel = @opWallInt.exactStokesDL;
    kernel2 = @opWallExt.exactStokesDL;
  end
  
  DLP = @(X) -1/2*X + opWallInt.exactStokesDLdiag(wallsInt,tt.wallDLPint,X);
  vback = opWallInt.nearSingInt(wallsInt,etaIntNew,DLP,[],NearWint2V,...
      kernel,kernelDirect,vesicle,false,false);
  
  DLP = @(X) -1/2*X + opWallExt.exactStokesDLdiag(wallsExt,tt.wallDLPext,X);
  vbackExt = opWallExt.nearSingInt(wallsExt,etaExtNew,DLP,[],NearWext2V,...
      kernel2,kernelDirect2,vesicle,false);
  vback = vback + vbackExt;
  
  for k = 2:nvbd
    stokeslet = RSnew(1:2,k); rotlet = RSnew(3,k);  
    vback = vback + tt.RSlets(vesicle.X,wallsInt.center(:,k-1),stokeslet,rotlet);
  end
end % if iconfined    

    
end % DNNsolveMixed2
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xnew,tenNew,etaIntNew,etaExtNew,RSnew] = DNNsolveMixed_wRigidDLP(o,Xold,...
        tenOld,etaIntOld,etaExtOld,RSold,Nnet,viscCont)

iuseTrueNear = true;

tt = o.tt;
op = tt.op;

opWallInt = tt.opWallInt;
opWallExt = tt.opWallExt;  
wallsInt = o.wallsInt;
wallsExt = o.wallsExt;
nvbdInt = wallsInt.nv;
nvbdExt = 1;
nvbd = nvbdExt + nvbdInt;
NbdInt = wallsInt.N;
NbdExt = wallsExt.N;
etaIntNew = zeros(2*NbdInt,nvbdInt);
etaExtNew = zeros(2*NbdExt,nvbdExt);
RSnew = zeros(3,nvbd);

% semi-implicit time stepping w/ splitting
vesicle = capsules(Xold,[],[],o.kappa,viscCont,1);
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

% Coefficient in the equation due to DLP
alpha = (viscCont(1)+1)/2;

% Build DLP only for the rigid particle
if viscCont(1)>1
  Nup = op.LPuprate * vesicle.N;
  Xup = [interpft(Xold(1:end/2,1),Nup);interpft(Xold(end/2+1:end,1),Nup)];
  vesicleUp = capsules(Xup,[],[],vesicle.kappa,viscCont(1),0); 
  DLPnoCorr = zeros(2*Nup,2*Nup,nv);
  DLPnoCorr(:,:,1) = op.stokesDLmatrixNoCorr(vesicleUp);
  
else
  DLPnoCorr = [];
end

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


% Compute bending forces + old tension forces
fBend = vesicle.tracJump(Xold,zeros(N,nv));
fTen = vesicle.tracJump(zeros(2*N,nv),tenOld);
tracJump = fBend+fTen;

% 1) EXPLICIT TENSION AT THE CURRENT STEP
if nv > 1
  
  if iuseTrueNear
  
  kernelDirect = @op.exactStokesSL;
  if tt.fmm
    kernel = @op.exactStokesSLfmm;
  else
    kernel = @op.exactStokesSL;
  end
  Galpert = op.stokesSLmatrix(vesicle);
  SLP = @(X) op.exactStokesSLdiag(vesicle,Galpert,X);
  farFieldtracJump = op.nearSingInt(vesicle,tracJump,SLP,SLPnoCorr,NearV2V,kernel,...
      kernelDirect,vesicle,true,false);
  
  else  
  % Far-field due to traction jump
  kernelDirect = @op.exactStokesSLregul;
  if tt.fmm
    kernel = @op.exactStokesSLfmm;
  else
    kernel = @op.exactStokesSLregul;
  end
  [~,nearField,farFieldtracJump] = op.divideNearFarSLP(vesicle,...
      tracJump,SLPnoCorr,NearV2V,kernel,kernelDirect,vesicle,true);      
  if o.useNear
    farFieldtracJump = farFieldtracJump + nearField;
  end
  end
else
  farFieldtracJump = zeros(2*N,nv);   
end

  

% Far-field due to walls or background flow
if o.confined
  if iuseTrueNear
  jump = -1/2;
  kernelDirect = @opWallInt.exactStokesDL;
  kernelDirect2 = @opWallExt.exactStokesDL;
  if tt.fmmDLP
    kernel = @opWallInt.exactStokesDLnewfmm;
    kernel2 = @opWallExt.exactStokesDLnewfmm;
  else
    kernel = @opWallInt.exactStokesDL;
    kernel2 = @opWallExt.exactStokesDL;
  end
  
  DLP = @(X) jump*X + opWallInt.exactStokesDLdiag(wallsInt,o.tt.wallDLPint,X);
  vback = opWallInt.nearSingInt(wallsInt,etaIntOld,DLP,[],NearWint2V,kernel,kernelDirect,vesicle,false,false);
  
  DLP = @(X) jump*X + opWallExt.exactStokesDLdiag(wallsExt,o.tt.wallDLPext,X);
  vbackExt = opWallExt.nearSingInt(wallsExt,etaExtOld,DLP,[],NearWext2V,kernel2,kernelDirect2,vesicle,false,false);
  vback = vback+vbackExt;
      
      
  else
  kernelDirect = @opWallInt.exactStokesDLregul;
  kernelDirect2 = @opWallExt.exactStokesDLregul;
  if tt.fmmDLP
    kernel = @opWallInt.exactStokesDLnewfmm;
    kernel2 = @opWallExt.exactStokesDLnewfmm;
  else
    kernel = @opWallInt.exactStokesDLregul;
    kernel2 = @opWallExt.exactStokesDLregul;
  end

  [~,nearField,vback] = opWallInt.divideNearFarSLP(wallsInt,etaIntOld,[],NearWint2V,...
      kernel,kernelDirect,vesicle,false);
  if o.useNear
    vback = vback + nearField;
  end
  
  [~,nearField,vbackExt] = opWallExt.divideNearFarSLP(wallsExt,etaExtOld,[],NearWext2V,...
      kernel2,kernelDirect2,vesicle,false);
  vback = vback + vbackExt;
  if o.useNear
    vback = vback + nearField;
  end
  end
  
  for k = 2:nvbd
    stokeslet = RSold(1:2,k); rotlet = RSold(3,k);  
    vback = vback + tt.RSlets(vesicle.X,wallsInt.center(:,k-1),stokeslet,rotlet);
  end
end



% if spherical particle has a visc. cont > 1, then first update it
if viscCont(1) > 1
  vesicleEll = capsules(Xold(:,1),[],[],vesicle.kappa,viscCont(1),1);
  vesicleEll.setUpRate();
  DLPEll = op.stokesDLmatrix(vesicleEll);
  DLPall = zeros(size(DLPEll,1),size(DLPEll,2),nv);
  DLPall(:,:,1) = DLPEll;
    
  % COMPUTE TENSION FOR RA = 0.9
  Nrigid = o.Nrigid;
  opRigid = o.opRigid;
  XrigUp = [interpft(Xold(1:end/2,1),Nrigid); interpft(Xold(end/2+1:end,1),Nrigid)];
  vesicleEllUp = capsules(XrigUp, [], [], vesicle.kappa, viscCont(1), 1);
  vesicleEllUp.setUpRate();
  SLPup = opRigid.stokesSLmatrix(vesicleEllUp);
  DLPup = opRigid.stokesDLmatrix(vesicleEllUp);
  DLPself = opRigid.exactStokesDLdiag(vesicleEllUp,DLPup,XrigUp);
  [BenUp,TenUp,DivUp] = vesicleEllUp.computeDerivs;

  vextEll = vback(:,1) + farFieldtracJump(:,1);
  vextEllUp = [interpft(vextEll(1:end/2),Nrigid);interpft(vextEll(end/2+1:end),Nrigid)];
  
  rhs = [vextEllUp * o.dt / alpha + XrigUp - DLPself / alpha; DivUp*XrigUp];
  matEll = [(eye(2*Nrigid) - DLPup/alpha) + ...
            o.dt/alpha*vesicle.kappa*SLPup*BenUp ...
        -o.dt/alpha*SLPup*TenUp; ...
        DivUp zeros(Nrigid)];
    
  solEll = matEll \ rhs;
  XnewEll = solEll(1:2*Nrigid);
  XnewEll = [interpft(XnewEll(1:end/2),N);interpft(XnewEll(end/2+1:end),N)];
  uEll = (XnewEll-Xold(:,1))/o.dt;
  tenNewEll = interpft(solEll(2*Nrigid+1:end),N);
  vesicleEll = capsules(XnewEll,[],[],vesicle.kappa,viscCont(1),1);
  vesicleEll.setUpRate();
  
  if iuseTrueNear
  % Velocity due to Viscosity Contrast (Old Config. is used)
  kernelDirect = @op.exactStokesDL;
  if tt.fmmDLP 
    kernel = @op.exactStokesDLnewfmm;
  else
    kernel = @op.exactStokesDL;
  end    
  jump = 1/2*(1-vesicle.viscCont);
  DLP = @(X) X*diag(jump) + op.exactStokesDLdiag(vesicle,DLPall,Xold);
  farFieldViscCont = op.nearSingInt(vesicle,[uEll zeros(2*N,nv-1)], DLP, DLPnoCorr, ...
      NearV2V,kernel, kernelDirect,vesicle,true,false);
  
  FDLPwallInt = op.nearSingInt(vesicle,[uEll zeros(2*N,nv-1)], DLP, [], NearV2Wint, ...
      kernel, kernelDirect, wallsInt, false, false);
  FDLPwallExt = op.nearSingInt(vesicle, [uEll zeros(2*N,nv-1)], DLP, [], NearV2Wext, ...
      kernel, kernelDirect, wallsExt, false, false);
  
  else
      
  % Velocity due to Viscosity Contrast (Old Config. is used)
  kernelDirect = @op.exactStokesDLregul;
  if tt.fmmDLP 
    kernel = @op.exactStokesDLnewfmm;
  else
    kernel = @op.exactStokesDLregul;
  end
  % On RBCs
  [~,nearField,farFieldViscCont] = op.divideNearFarSLP(vesicle,...
      [Xold(:,1) zeros(2*N,nv-1)],DLPnoCorr,NearV2V,kernel,kernelDirect,vesicle,true);
  if o.useNear
    farFieldViscCont = farFieldViscCont + nearField; 
  end
  
  % On walls
  [~,nearField,FDLPwallInt] = op.divideNearFarSLP(vesicle,[Xold(:,1) zeros(2*N,nv-1)],DLPnoCorr,...
      NearV2Wint,kernel,kernelDirect,wallsInt,false);
  if o.useNear
    FDLPwallInt = FDLPwallInt + nearField;
  end
  
  [~,nearField,FDLPwallExt] = op.divideNearFarSLP(vesicle,[Xold(:,1) zeros(2*N,nv-1)],DLPnoCorr,...
      NearV2Wext,kernel,kernelDirect,wallsExt,false);
  if o.useNear
    FDLPwallExt = FDLPwallExt + nearField;
  end
  end
  
else
  farFieldViscCont = zeros(2*N,nv);
  FDLPwallInt = zeros(2*NbdInt,nvbdInt);
  FDLPwallExt = zeros(2*NbdExt,nvbdExt);
end



scaling = zeros(nv-1,1); rotate = zeros(nv-1,1); trans = zeros(2,nv-1);
sortIdx = zeros(Nnet,nv-1);
XinputFourier = zeros(o.nComp,1,1,nv-1);
XinputTenBen = zeros(o.nComp,1,1,nv-1);
for k = 2 : nv
  % Standardize vesicle  
  [Xstand,scaling(k-1),rotate(k-1),trans(:,k-1),sortIdx(:,k-1)] = ...
      o.standardizationStep(Xold(:,k),Nnet);
  
  % Prepare inputs for Fourier network and self-bending network
  XinputFourier(:,:,:,k-1) = o.prepareInputForNet(Xstand,'tensionOnFourier');
  XinputTenBen(:,:,:,k-1) = o.prepareInputForNet(Xstand,'tensionSelfBend');

end
% Approximate inv(Div*G*Ten)*Div*vExt
vBackSolve = o.invTenMatOnVback(XinputFourier,vback(:,2:end) + farFieldtracJump(:,2:end) + farFieldViscCont(:,2:end),...
    scaling,rotate,sortIdx);
% Approximate inv(Div*G*Ten)*G*(-Ben)*x
selfBendSolve = o.invTenMatOnSelfBend(XinputTenBen,scaling,sortIdx,N);
% update tension
tenNew = -(vBackSolve+selfBendSolve);

if viscCont(1) > 1
  tenNew = [tenNewEll tenNew];
  tracJumpEll = vesicleEll.tracJump(XnewEll,tenNewEll);    
  % Update the traction jump
  fTen = vesicle.tracJump(zeros(2*N,nv),tenNew); tracJump = fBend+fTen;
  tracJump(:,1) = tracJumpEll;
else
  % Update the traction jump
  fTen = vesicle.tracJump(zeros(2*N,nv),tenNew); tracJump = fBend+fTen;
end



if o.confined
  % 2) SOLVE FOR DENSITY and RS ON WALLS
  % vesicle2wall interactions
  
  if iuseTrueNear
  kernelDirect = @op.exactStokesSL;
  if tt.fmm
    kernel = @op.exactStokesSLfmm;
  else
    kernel = @op.exactStokesSL;
  end
  
  SLP = @(X) op.exactStokesSLdiag(vesicle,Galpert,X);
  ves2wallInt = op.nearSingInt(vesicle,tracJump,SLP,[],NearV2Wint,kernel,...
      kernelDirect,wallsInt,false,false);
  ves2wallExt = op.nearSingInt(vesicle,tracJump,SLP,[],NearV2Wext,kernel,...
      kernelDirect,wallsExt,false,false);
      
  else
  kernelDirect = @op.exactStokesSLregul;
  if tt.fmm
    kernel = @op.exactStokesSLfmm;
  else
    kernel = @op.exactStokesSLregul;
  end
  
  [~,nearField,ves2wallInt] = op.divideNearFarSLP(vesicle,tracJump,[],NearV2Wint,...
      kernel,kernelDirect,wallsInt,0);
  if o.useNear
    ves2wallInt = ves2wallInt+nearField;
  end
 
  [~,nearField,ves2wallExt] = op.divideNearFarSLP(vesicle,tracJump,[],NearV2Wext,...
      kernel,kernelDirect,wallsExt,0);
  if o.useNear
    ves2wallExt = ves2wallExt+nearField;
  end
  end
  
  RHSint = wallsInt.u(:)-ves2wallInt(:)-FDLPwallInt(:);
  RHSext = wallsExt.u(:)-ves2wallExt(:)-FDLPwallExt(:);
  RHS = [RHSext;RHSint;zeros(3*(nvbd-1),1)];
  etaRS = tt.bdiagWall*RHS;
  etaExtNew = etaRS(1:2*NbdExt);
  etaRS = etaRS(2*NbdExt+1:end);
  for iw = 1 : nvbdInt
    etaIntNew(:,iw) = etaRS((iw-1)*2*NbdInt+1:iw*2*NbdInt);  
    if iw <= nvbdInt
      RSnew(:,iw+1) = etaRS(2*NbdInt*nvbdInt+(iw-1)*3+1:2*NbdInt*nvbdInt+iw*3);
    end
  end 
  
  
  
  if iuseTrueNear
  jump = -1/2;    
  % Update vback due to new density,rotlets and stokeslet
  kernelDirect = @opWallInt.exactStokesDL;
  kernelDirect2 = @opWallExt.exactStokesDL;
  if tt.fmmDLP
    kernel = @opWallInt.exactStokesDLnewfmm;
    kernel2 = @opWallExt.exactStokesDLnewfmm;
  else
    kernel = @opWallInt.exactStokesDL;
    kernel2 = @opWallExt.exactStokesDL;
  end    
  DLP = @(X) jump*X + opWallInt.exactStokesDLdiag(wallsInt,tt.wallDLPint,X);
  vback = opWallInt.nearSingInt(wallsInt,etaIntNew,DLP,[],NearWint2V,kernel,kernelDirect,vesicle,false,false);
  
  DLP = @(X) jump*X + opWallExt.exactStokesDLdiag(wallsExt,tt.wallDLPext,X);
  vbackExt = opWallExt.nearSingInt(wallsExt,etaExtNew,DLP,[],NearWext2V,kernel2,kernelDirect2,vesicle,false,false);
  
  vback = vback + vbackExt;
      
  else
  % Update vback due to new density,rotlets and stokeslet
  kernelDirect = @opWallInt.exactStokesDLregul;
  kernelDirect2 = @opWallExt.exactStokesDLregul;
  if tt.fmmDLP
    kernel = @opWallInt.exactStokesDLnewfmm;
    kernel2 = @opWallExt.exactStokesDLnewfmm;
  else
    kernel = @opWallInt.exactStokesDLregul;
    kernel2 = @opWallExt.exactStokesDLregul;
  end
  
 
  [~,nearField,vback] = opWallInt.divideNearFarSLP(wallsInt,etaIntNew,[],NearWint2V,...
      kernel,kernelDirect,vesicle,false);
  if o.useNear
    vback = vback + nearField;
  end
 
  [~,nearField,vbackExt] = opWallExt.divideNearFarSLP(wallsExt,etaExtNew,[],NearWext2V,...
      kernel2,kernelDirect2,vesicle,false);
  vback = vback + vbackExt;
  if o.useNear
    vback = vback + nearField;
  end
  end
  
  for k = 2:nvbd
    stokeslet = RSnew(1:2,k); rotlet = RSnew(3,k);  
    vback = vback + tt.RSlets(vesicle.X,wallsInt.center(:,k-1),stokeslet,rotlet);
  end
end % if iconfined    

% 3) SOLVE FOR NEXT TIME STEP USING OPERATOR SPLITTING
% EXPLICIT-TREATMENT OF INTER-VESICLE BENDING AND TENSION
% TENSION AND X are coupled (semi-implicit).
if nv > 1
  if iuseTrueNear
  kernelDirect = @op.exactStokesSL;
  if tt.fmm
    kernel = @op.exactStokesSLfmm;
  else
    kernel = @op.exactStokesSL;
  end
  SLP = @(X) op.exactStokesSLdiag(vesicle,Galpert,X);
  farFieldtracJump = op.nearSingInt(vesicle,tracJump,SLP,SLPnoCorr,NearV2V,kernel,...
      kernelDirect,vesicle,true,false);
  else
  % compute traction jump due to explicit tension and old bending
  kernelDirect = @op.exactStokesSLregul;
  if tt.fmm
    kernel = @op.exactStokesSLfmm;
  else
    kernel = @op.exactStokesSLregul;
  end

  [~,nearField,farFieldtracJump] = op.divideNearFarSLP(...
          vesicle,tracJump,SLPnoCorr,NearV2V,kernel,kernelDirect,vesicle,1);
  if o.useNear
    farFieldtracJump = farFieldtracJump + nearField;
  end    
  end
else
  farFieldtracJump = zeros(2*N,nv);  
end

% Walls2Ves is already computed if confined (vback)    
rotate = zeros(nv-1,1); trans = zeros(2,nv-1); sortIdx = zeros(Nnet,nv-1);
scaling = zeros(nv-1,1); Xinput = zeros(o.nComp,1,1,nv-1);
for k = 2 : nv
  % 1) TRANSLATION W/ NETWORK
  % Standardize vesicle
  [Xstand,scaling(k-1),rotate(k-1),trans(:,k-1),sortIdx(:,k-1)] = ...
      o.standardizationStep(Xold(:,k),Nnet);
  % Prepare input for advection network
  Xinput(:,1,1,k-1) = o.prepareInputForNet(Xstand,'advection');
end
% Take a step due to advaction
Xmid = o.translateVinfwNN(Xinput,vback(:,2:end) + farFieldtracJump(:,2:end) + farFieldViscCont(:,2:end),...
  Xold(:,2:end),rotate,sortIdx);


rotate = zeros(nv-1,1); trans = zeros(2,nv-1); sortIdx = zeros(Nnet,nv-1);
scaling = zeros(nv-1,1); Xinput = zeros(o.nCompRelax,1,1,nv-1);
for k = 2 : nv
  % 2) RELAXATION w/ NETWORK
  % Standardize vesicle Xmid
  [Xstand,scaling(k-1),rotate(k-1),trans(:,k-1),sortIdx(:,k-1)] = ...
      o.standardizationStep(Xmid(:,k-1),Nnet);
  % Prepare input for relaxation network
  Xinput(:,1,1,k-1) = o.prepareInputForNet(Xstand,'relaxation');
end
% Take a step
Xnew = o.relaxWNN(Xinput,scaling,rotate,trans,sortIdx,N);

% Add rigid particle
Xnew = [XnewEll Xnew];
    
end % DNNsolveMixed_wRigidDLP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xnew,tenNew,etaIntNew,etaExtNew,RSnew] = DNNsolveMixed_wRigidDLP2(o,Xold,...
        tenOld,etaIntOld,etaExtOld,RSold,Nnet,viscCont)
% USING TRUE NEAR INTERACTIONS

tt = o.tt;
op = tt.op;
opRigid = o.opRigid;
Nrigid = o.Nrigid;

opWallInt = tt.opWallInt;
opWallExt = tt.opWallExt;  
wallsInt = o.wallsInt;
wallsExt = o.wallsExt;
nvbdInt = wallsInt.nv;
nvbdExt = 1;
nvbd = nvbdExt + nvbdInt;
NbdInt = wallsInt.N;
NbdExt = wallsExt.N;
etaIntNew = zeros(2*NbdInt,nvbdInt);
RSnew = zeros(3,nvbd);

% Vesicle class with old configuration
vesicle = capsules(Xold,[],[],o.kappa,viscCont,1);
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

% Coefficient in the equation due to DLP for Rigid only
alpha = (viscCont(1)+1)/2;

% Build DLP only for the rigid particle
if viscCont(1)>1
  Nup = op.LPuprate * vesicle.N;
  Xup = [interpft(Xold(1:end/2,1),Nup);interpft(Xold(end/2+1:end,1),Nup)];
  vesicleUp = capsules(Xup,[],[],vesicle.kappa,viscCont(1),0); 
  DLPnoCorr = zeros(2*Nup,2*Nup,nv);
  DLPnoCorr(:,:,1) = op.stokesDLmatrixNoCorr(vesicleUp);
else
  DLPnoCorr = [];
end

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

% Compute bending forces + old tension forces
fBend = vesicle.tracJump(Xold,zeros(N,nv));
fTen = vesicle.tracJump(zeros(2*N,nv),tenOld);
tracJump = fBend+fTen;

% 1) EXPLICITLY SOLVE FOR Xn+1

% far-field due to traction jump
kernelDirect = @op.exactStokesSL;
if tt.fmm
  kernel = @op.exactStokesSLfmm;
else
  kernel = @op.exactStokesSL;
end
Galpert = op.stokesSLmatrix(vesicle);
SLP = @(X) op.exactStokesSLdiag(vesicle,Galpert,X);
farFieldtracJump = op.nearSingInt(vesicle,tracJump,SLP,SLPnoCorr,NearV2V,kernel,...
  kernelDirect,vesicle,true,false);

% far-field due to walls
jump = -1/2;
kernelDirect = @opWallInt.exactStokesDL;
kernelDirect2 = @opWallExt.exactStokesDL;
if tt.fmmDLP
  kernel = @opWallInt.exactStokesDLnewfmm;
  kernel2 = @opWallExt.exactStokesDLnewfmm;
else
  kernel = @opWallInt.exactStokesDL;
  kernel2 = @opWallExt.exactStokesDL;
end

DLP = @(X) jump*X + opWallInt.exactStokesDLdiag(wallsInt,o.tt.wallDLPint,X);
vback = opWallInt.nearSingInt(wallsInt,etaIntOld,DLP,[],NearWint2V,kernel,kernelDirect,vesicle,false,false);

DLP = @(X) jump*X + opWallExt.exactStokesDLdiag(wallsExt,o.tt.wallDLPext,X);
vbackExt = opWallExt.nearSingInt(wallsExt,etaExtOld,DLP,[],NearWext2V,kernel2,kernelDirect2,vesicle,false,false);
vback = vback+vbackExt;

for k = 2:nvbd
  stokeslet = RSold(1:2,k); rotlet = RSold(3,k);  
  vback = vback + tt.RSlets(vesicle.X,wallsInt.center(:,k-1),stokeslet,rotlet);
end

% far-field due to viscosity contrast
if viscCont(1) > 1
  vesicleEll = capsules(Xold(:,1),[],[],vesicle.kappa,viscCont(1),1);
  vesicleEll.setUpRate();
  DLPEll = op.stokesDLmatrix(vesicleEll);
  DLPall = zeros(size(DLPEll,1),size(DLPEll,2),nv);
  DLPall(:,:,1) = DLPEll;
    
  XrigUp = [interpft(Xold(1:end/2,1),Nrigid); interpft(Xold(end/2+1:end,1),Nrigid)];
  vesicleEllUp = capsules(XrigUp, [], [], vesicle.kappa, viscCont(1), 1);
  vesicleEllUp.setUpRate();
  SLPup = opRigid.stokesSLmatrix(vesicleEllUp);
  DLPup = opRigid.stokesDLmatrix(vesicleEllUp);
  DLPself = opRigid.exactStokesDLdiag(vesicleEllUp,DLPup,XrigUp);
  [BenUp,TenUp,DivUp] = vesicleEllUp.computeDerivs;

  vextEll = vback(:,1) + farFieldtracJump(:,1);
  vextEllUp = [interpft(vextEll(1:end/2),Nrigid);interpft(vextEll(end/2+1:end),Nrigid)];
  
  % Solve for rigid body
  if 1
  rhs = [vextEllUp * o.dt / alpha + XrigUp - DLPself / alpha; DivUp*XrigUp];

  matEll = [(eye(2*Nrigid) - DLPup/alpha) + ...
            o.dt/alpha*vesicle.kappa*SLPup*BenUp ...
        -o.dt/alpha*SLPup*TenUp; ...
        DivUp zeros(Nrigid)];
    
  else
  solveIminusD = (alpha*eye(2*Nrigid) - DLPup) \ vextEllUp;

  rhs = [vextEllUp * o.dt / alpha + XrigUp - DLPself / alpha; -DivUp*(solveIminusD)];
  
  matEll = [(eye(2*Nrigid) - DLPup/alpha) + ...
            o.dt/alpha*vesicle.kappa*SLPup*BenUp ...
            -o.dt/alpha*SLPup*TenUp; ...
            -vesicle.kappa*DivUp*((alpha*eye(2*Nrigid) - DLPup)\(SLPup*BenUp)) ...
            DivUp*((alpha*eye(2*Nrigid) - DLPup)\eye(2*Nrigid))*(SLPup*TenUp)];
  end
  solEll = matEll \ rhs;
  XnewEll = solEll(1:2*Nrigid);
  XnewEll = [interpft(XnewEll(1:end/2),N);interpft(XnewEll(end/2+1:end),N)];
  uRigid = (XnewEll-Xold(:,1))/o.dt;
  tenNewEll = interpft(solEll(2*Nrigid+1:end),N);
  vesicleEll = capsules(XnewEll,[],[],vesicle.kappa,viscCont(1),1);
  vesicleEll.setUpRate();    
  
  % Velocity due to Viscosity Contrast (Old Config. is used)
  kernelDirect = @op.exactStokesDL;
  if tt.fmmDLP 
    kernel = @op.exactStokesDLnewfmm;
  else
    kernel = @op.exactStokesDL;
  end    
  jump = 1/2*(1-vesicle.viscCont);
  DLP = @(X) X*diag(jump) + op.exactStokesDLdiag(vesicle,DLPall,X);
  farFieldViscCont = op.nearSingInt(vesicle,[uRigid zeros(2*N,nv-1)], DLP, DLPnoCorr, ...
      NearV2V,kernel, kernelDirect,vesicle,true,false);
  
  FDLPwallInt = op.nearSingInt(vesicle, [uRigid zeros(2*N,nv-1)], DLP, [], NearV2Wint, ...
      kernel, kernelDirect, wallsInt, false, false);
  FDLPwallExt = op.nearSingInt(vesicle, [uRigid zeros(2*N,nv-1)], DLP, [], NearV2Wext, ...
      kernel, kernelDirect, wallsExt, false, false);

end


% Walls2Ves is already computed if confined (vback)    
rotate = zeros(nv-1,1); trans = zeros(2,nv-1); sortIdx = zeros(Nnet,nv-1);
scaling = zeros(nv-1,1); Xinput = zeros(o.nComp,1,1,nv-1);
for k = 2 : nv
  % 1) TRANSLATION W/ NETWORK
  % Standardize vesicle
  [Xstand,scaling(k-1),rotate(k-1),trans(:,k-1),sortIdx(:,k-1)] = ...
      o.standardizationStep(Xold(:,k),Nnet);
  % Prepare input for advection network
  Xinput(:,1,1,k-1) = o.prepareInputForNet(Xstand,'advection');
end
% Take a step due to advaction
Xmid = o.translateVinfwNN(Xinput,vback(:,2:end) + farFieldtracJump(:,2:end) + farFieldViscCont(:,2:end),...
  Xold(:,2:end),rotate,sortIdx);


rotate = zeros(nv-1,1); trans = zeros(2,nv-1); sortIdx = zeros(Nnet,nv-1);
scaling = zeros(nv-1,1); Xinput = zeros(o.nCompRelax,1,1,nv-1);
for k = 2 : nv
  % 2) RELAXATION w/ NETWORK
  % Standardize vesicle Xmid
  [Xstand,scaling(k-1),rotate(k-1),trans(:,k-1),sortIdx(:,k-1)] = ...
      o.standardizationStep(Xmid(:,k-1),Nnet);
  % Prepare input for relaxation network
  Xinput(:,1,1,k-1) = o.prepareInputForNet(Xstand,'relaxation');
end
% Take a step
Xnew = o.relaxWNN(Xinput,scaling,rotate,trans,sortIdx,N);

% Add rigid particle
Xnew = [XnewEll Xnew];

% 2) COMPUTE TENSION NOW
% Update bending SLP
vesicleNew = capsules(Xnew,[],[],vesicle.kappa,viscCont,1);
vesicleNew.setUpRate();
fBend = vesicleNew.tracJump(Xnew,zeros(N,nv));

farFieldBending = op.nearSingInt(vesicle,fBend,SLP,SLPnoCorr,NearV2V,kernel,...
      kernelDirect,vesicle,true,false);
  
scaling = zeros(nv-1,1); rotate = zeros(nv-1,1); trans = zeros(2,nv-1);
sortIdx = zeros(Nnet,nv-1);
XinputFourier = zeros(o.nComp,1,1,nv-1);
XinputTenBen = zeros(o.nComp,1,1,nv-1);
for k = 2 : nv
  % Standardize vesicle  
  [Xstand,scaling(k-1),rotate(k-1),trans(:,k-1),sortIdx(:,k-1)] = ...
      o.standardizationStep(Xnew(:,k),Nnet);
  
  % Prepare inputs for Fourier network and self-bending network
  XinputFourier(:,:,:,k-1) = o.prepareInputForNet(Xstand,'tensionOnFourier');
  XinputTenBen(:,:,:,k-1) = o.prepareInputForNet(Xstand,'tensionSelfBend');

end
% Approximate inv(Div*G*Ten)*Div*vExt
vBackSolve = o.invTenMatOnVback(XinputFourier,vback(:,2:end) + farFieldBending(:,2:end) + farFieldViscCont(:,2:end),...
    scaling,rotate,sortIdx);
% Approximate inv(Div*G*Ten)*G*(-Ben)*x
selfBendSolve = o.invTenMatOnSelfBend(XinputTenBen,scaling,sortIdx,N);
% update tension
tenNew = -(vBackSolve+selfBendSolve);
tenNew = [tenNewEll tenNew];

% Update the traction jump
fTen = vesicleNew.tracJump(zeros(2*N,nv),tenNew); tracJump = fBend+fTen;

  
% Compute velocity on the walls based on the newest X and tension
kernelDirect = @op.exactStokesSL;
if tt.fmm
  kernel = @op.exactStokesSLfmm;
else
  kernel = @op.exactStokesSL;
end
  
SLP = @(X) op.exactStokesSLdiag(vesicle,Galpert,X);
ves2wallInt = op.nearSingInt(vesicle,tracJump,SLP,[],NearV2Wint,kernel,...
      kernelDirect,wallsInt,false,false);
ves2wallExt = op.nearSingInt(vesicle,tracJump,SLP,[],NearV2Wext,kernel,...
      kernelDirect,wallsExt,false,false);

RHSint = wallsInt.u(:)-ves2wallInt(:)-FDLPwallInt(:);
RHSext = wallsExt.u(:)-ves2wallExt(:)-FDLPwallExt(:);
RHS = [RHSext;RHSint;zeros(3*(nvbd-1),1)];
etaRS = tt.bdiagWall*RHS;
etaExtNew = etaRS(1:2*NbdExt);
etaRS = etaRS(2*NbdExt+1:end);
for iw = 1 : nvbdInt
  etaIntNew(:,iw) = etaRS((iw-1)*2*NbdInt+1:iw*2*NbdInt);  
  if iw <= nvbdInt
    RSnew(:,iw+1) = etaRS(2*NbdInt*nvbdInt+(iw-1)*3+1:2*NbdInt*nvbdInt+iw*3);
  end
end 
end % DLP2
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
prams.totnv = pramsIn.totnv;
prams.Nbd = 0;
prams.N = N;
prams.nv = nv;
prams.NbdInt = pramsIn.NbdInt;
prams.NbdExt = pramsIn.NbdExt;
prams.nvbd = pramsIn.nvbdInt+pramsIn.nvbdExt;
prams.nvbdInt = pramsIn.nvbdInt; prams.nvbdExt = pramsIn.nvbdExt;
prams.kappa = pramsIn.kappa;
prams.nrow = pramsIn.nrow;
prams.ncol = pramsIn.ncol;
prams.Dpostx = pramsIn.Dpostx;
prams.Dposty = pramsIn.Dposty;
prams.Dx = pramsIn.Dx; 
prams.Dy = pramsIn.Dy;
prams.epsilon = pramsIn.epsilon;

options.diffDiscWalls = 1;

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
options.repulsion = o.repulsion;

options.confined = true;
options.farField = pramsIn.farField;
options.farFieldSpeed = pramsIn.speed;

[options,prams] = initVes2D(options,prams);
    
om = monitor(X,options,prams);
tt = tstep(options,prams,om);

if strcmp(pramsIn.farField,'rotDLD')
  tt.farField = @(X,Xint) tt.bgFlow(X,options.farField,...
      'intWalls',Xint,'velG',pramsIn.gPer,'Speed',1);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % methods

end % dnnTools
