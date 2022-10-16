clear; clc;
addpath ../src/
addpath ../examples/
disp('Single vesicle in DLD')

% DLD parameters
%-------------------------------------------------------------------------
pillarCross = 'circlePillar'; % or trianglePillar
Dx = 12; % gap size in x-dir. (micron)
Dy = 12; % gap size in y-dir. (micron)
runName = 'ra90_kb10_circPost';
wbox = 25; % square box size around a pillar (including Dx and Dy) (micron)
theta = atan(1/6); 
prams.NbdExt = 3072; % # of points on exterior wall % 3072 for 6 rows
prams.NbdInt = 32; % # of points on posts % 64
prams.N = 32; % # of points per vesicle
prams.speed = 100; % Maximum velocity in the gap scale
prams.vesViscCont = 1; % viscosity contrast
prams.kappa = 10; % bending rigidity
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
exactSolveFreq = 1; % solve exactly at every [] time steps
prams.farField = 'DLD';
prams.Th = 1.5; % time horizon
prams.nv = 1; prams.totnv = prams.nv;
prams.fmm = useFMM; % use FMM for ves2ves
prams.fmmDLP = useFMM; % use FMM for ves2walls
prams.dt = 1E-4; % time step size

% NO NEED TO CHANGE ANYTHING BELOW THE LINE
%-------------------------------------------------------------------------
oc = curve;
% WALLS and VESICLES:
% -------------------------------------------------------------------------
% DLD Geometry:
scale = 0.05;
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
prams.ncol = ceil(1.5*prams.periodN);
% number of obstacles in y direction (number of rows)
prams.nrow    = 6; 

[XwallsInt,XwallsExt,prams.Lext,prams.xfreeze,radiusRight] = ...
    initializeDLD(prams,Xpost,pillarCross);
prams.nvbdInt = size(XwallsInt,2);
prams.nvbdExt = 1;

% Vesicle:
if 0
load relaxedVesIC
X = [interpft(X(1:end/2),prams.N);interpft(X(end/2+1:end),prams.N)];
end
X = oc.ellipse(prams.N,0.9); [~,a1,l1] = oc.geomProp(X); X = X/l1;

IA = oc.getIncAngle(X);
% start at the vertical location parallel to the 2nd pillar
% centy = -prams.Dy*0.45; /
centy = -prams.Dy*0.2;
centx = min(XwallsExt(1:end/2))+radiusRight+0.6*prams.Dy; 
angles = 0*ones(prams.nv,1);  

Xnew = zeros(size(X));
Xnew(1:end/2) = cos(-IA+angles)*X(1:end/2)-sin(-IA+angles)*X(end/2+1:end);
Xnew(end/2+1:end) = sin(-IA+angles)*X(1:end/2)+cos(-IA+angles)*X(end/2+1:end);
X = [Xnew(1:end/2)-mean(Xnew(1:end/2))+centx;...
      Xnew(end/2+1:end)-mean(Xnew(end/2+1:end))+centy];
[~,area0,len0] = oc.geomProp(X);
% viscosity contrast
prams.viscCont = prams.vesViscCont*ones(prams.nv,1);  

if 0
figure(1);clf;hold on;
plot([XwallsExt(1:end/2);XwallsExt(1)],...
    [XwallsExt(end/2+1:end);XwallsExt(end/2+1)],'k','linewidth',3)
plot([XwallsInt(1:end/2,:);XwallsInt(1,:)],[XwallsInt(end/2+1:end,:);...
    XwallsInt(end/2+1,:)],'k','linewidth',2)
plot([X(1:end/2);X(1)],[X(end/2+1:end);X(end/2+1)],'r','linewidth',2)
axis equal
pause
end
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
%-------------------------------------------------------------------------
disp(['N = ' num2str(prams.N) ', nv = ' num2str(prams.nv) ', Th = ' num2str(prams.Th)])
%-------------------------------------------------------------------------

% BUILD DNN CLASS
% -------------------------------------------------------------------------
dnn = dnnToolsDLD(X,XwallsInt,XwallsExt,prams);
tt = dnn.tt; dnn.oc = oc; dnn.useNear = iuseNear; dnn.useTrueNear = iTrueNear;
dnn.repulsion = irepulsion;
if irepulsion
  vesicle = capsules(X(:,1),[],[],prams.kappa,1,true);
  vesicle.setUpRate();
  dnn.minDist = 0.4;
  dnn.repStrength = vesicle.repulStrengthScale(dnn.minDist,tt,prams.speed);
end    
% -------------------------------------------------------------------------

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
Xhist = zeros(2*prams.N,prams.nv,ntime); Xhist(:,:,1) = X; 
sigStore = zeros(prams.N,prams.nv,ntime); sigStore(:,:,1) = sig;
errALPred = zeros(ntime,1);

ncountCNN = 0;
ncountExct = 0;
% ------------------------------------------------------------------------
fileName = ['./output/' runName 'DNN'];


% TIME STEPPING
for it = 2 : ntime
  disp('********************************************') 
  disp([num2str(it) 'th (/' num2str(ntime) ...
    ') time step, time: ' num2str(time(it))])

  
  if rem(it,exactSolveFreq) == 0      
    % EXACT SOLVE  
    disp('Taking an exact time step...'); tStart = tic;  
    ncountExct = ncountExct + 1;
    vesicle = capsules(Xhist(:,:,it-1),[],[],prams.kappa,ones(prams.nv,1),1);
    vesicle.setUpRate();
    [Xhist(:,:,it),sigStore(:,:,it),etaInt,etaExt,RS,...
        iter,iflag] = tt.timeStepSimpleDiffDisc(Xhist(:,:,it-1),...
        sigStore(:,:,it-1),etaInt,etaExt,RS,ones(prams.nv,1),...
        dnn.wallsInt,dnn.wallsExt,vesicle);
    
  else 
    % NETWORK SOLVE 
    ncountCNN = ncountCNN+1;
    disp('Taking a step with DNNs...');  tStart = tic;
    [Xhist(:,:,it),sigStore(:,:,it),etaInt,etaExt,RS] = dnn.DNNsolve(...
        Xhist(:,:,it-1),sigStore(:,:,it-1),etaInt,etaExt,RS,Nnet);
    
    % check if there is any self-intersecting shape
    for k = 1 : prams.nv
      [xIntersect,~,~] = oc.selfintersect(Xhist(:,k,it));
      if ~isempty(xIntersect)
        disp('New vesicle shape is self-intersecting!!!')
      end
    end

    if ~isempty(xIntersect) 
      disp('Taking an exact time step...')
      ncountCNN = ncountCNN-1;
      ncountExct = ncountExct+1;
      vesicle = capsules(Xhist(:,:,it-1),[],[],prams.kappa,ones(prams.nv,1),1);
      vesicle.setUpRate();
      [Xhist(:,:,it),sigStore(:,:,it),etaInt,etaExt,RS,...
        iter,iflag] = tt.timeStepSimpleDiffDisc(Xhist(:,:,it-1),...
        sigStore(:,:,it-1),etaInt,etaExt,RS,ones(prams.nv,1),...
        dnn.wallsInt,dnn.wallsExt,vesicle);
    end % self-intersection occurs
    
  end % if rem(it,10) == 0

  % JIGGLING
  if iJiggle % jiggle vesicles pointwise if there is a near-collision
    Xhist(:,:,it) = oc.fixCollisionsDLD(Xhist(:,:,it),XwallsInt);

    % Equally distribute points in arc-length b/c jiggling messes up points
    disp('Equally distributing points on arc-length...') 
    Xiter = Xhist(:,:,it);
    for iter = 1 : 5
      [Xiter,~,~] = oc.redistributeArcLength(Xiter);
    end

    % Fix misalignment in center and angle due to reparametrization
    Xhist(:,:,it) = oc.alignCenterAngle(Xhist(:,:,it),Xiter);

  end
  
  % AREA-LENGTH CORRECTION
  disp('Correcting area-length errors...')
  [Xnew,ifail] = oc.correctAreaAndLength2(Xhist(:,:,it),area0,len0);
  
  if ifail
    disp('Error in AL cannot be corrected!!!')
  end
  Xhist(:,:,it) = oc.alignCenterAngle(Xhist(:,:,it),Xnew);

  % check if shape is intersecting
  for k = 1 : prams.nv
    [xIntersect,~,~] = oc.selfintersect(Xhist(:,k,it));
    if ~isempty(xIntersect)
      disp('New vesicle shape is self-intersecting!!!')
    end
  end
  
  % Equally distribute points in arc-length
  disp('Equally distributing points on arc-length...')
  Xiter = Xhist(:,:,it);
  for iter = 1 : 5
    [Xiter,~,~] = oc.redistributeArcLength(Xiter);
  end

  % Fix misalignment in center and angle due to reparametrization
  Xhist(:,:,it) = oc.alignCenterAngle(Xhist(:,:,it),Xiter);
  
  disp(['took ' num2str(toc(tStart)) ' seconds.'])
  
  % Compute error in area and length
  [~,area,len] = oc.geomProp(Xhist(:,:,it));
  errArea = max(abs(area-area0)./area0); errLen = max(abs(len-len0)./len0);
  errALPred(it) = max(errArea,errLen);
  disp(['Error in area and length: ' num2str(errALPred(it))])   
  disp('********************************************') 
  disp(' ')
  
  % Check whether vesicle should be frozen
  if mean(Xhist(1:end/2,1,it)) >= prams.xfreeze
    disp('Freeze vesicle and terminate simulation...')
    save(fileName,'Xhist','time','sigStore','it',...
       'errALPred','ncountCNN','ncountExct','XwallsInt','XwallsExt')
   break
  end % if mean
end

% Save data to a mat-file:
save(fileName,'Xhist','time','sigStore','it',...
'errALPred','ncountCNN','ncountExct','XwallsInt','XwallsExt')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [XwallsInt,XwallsExt,L,xfreeze,radiusRight] = ...
    initializeDLD(prams,Xpost,pillarCross)
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

% Height and length of the exterior wall
H = (prams.nrow+1)*prams.Dy + prams.nrow*prams.Dposty; % necessary for large periodN
L = (prams.ncol+1)*(prams.Dx+prams.Dpostx);

% parameters for the exterior wall
a = L/2;
b = H/2; order = 20;


% a and b control length and height of exterior wall.
% Order controls the regularity.
t = (0:prams.NbdExt-1)'*2*pi/prams.NbdExt;
r = (cos(t).^order + sin(t).^order).^(-1/order);
x = a*r.*cos(t); y = b*r.*sin(t);
XwallsExt = [x;y];

% Row shift  
delLat = (prams.Dy+prams.Dposty)*prams.epsilon;

% Generate more rows than asked and remove ones outside exterior wall
pnrow = ceil(4*prams.nrow);
pHeight = (pnrow-1)*prams.Dy + (pnrow-1)*prams.Dposty;

% Place horizontally the obstacles in the first column
centy1stCol = linspace(-pHeight/2,pHeight/2,pnrow)';
centx = linspace(-L/2+1.5*prams.Dpostx-maxDist2Left+prams.Dx,...
    L/2-prams.Dpostx/2-maxDist2Left-prams.Dx,prams.ncol);

% when to freeze vesicles (after slightly after one period)
if numel(centx)>=ceil(prams.periodN)+2
  xfreeze = centx(1*ceil(prams.periodN)+2)+prams.Dpostx/2;
else
  xfreeze = [];
end

% Place vertically all the obstacles
centx = centx(ones(pnrow,1),:);

% Shift centers of the obstacles
delLatVect = [0:prams.ncol-1]*delLat;
centy = delLatVect(ones(pnrow,1),:) + centy1stCol(:,ones(1,prams.ncol));

% Place the obstacles
Xobs = zeros(2*prams.NbdInt,pnrow*prams.ncol);
for iwall = 1 : pnrow*prams.ncol
  Xobs(:,iwall) = [Xpost(1:end/2)+centx(iwall);...
    Xpost(end/2+1:end)+centy(iwall)];
end


% Remove obstacles outside the geometry, replace the ones intersecting the
% geometry with ellipses
tobs = [0:prams.NbdInt-1]'*2*pi/prams.NbdInt;
jdx = 1;
for iwall = 1 : pnrow*prams.ncol
  xwall = Xobs(1:end/2,iwall); ywall = Xobs(end/2+1:end,iwall);
  idx1 = isempty(find(abs(xwall) >= 0.98*L,1));
  idx2 = isempty(find(abs(ywall) >= 0.98*H/2,1));
  if idx1 && idx2
    XwallsInt(:,jdx) = [xwall;ywall];
    jdx = jdx + 1;
  else
    centx = mean(xwall);
    if isempty(find(ywall<0,1))
      minyloc = min(ywall);
      smallDiam = 0.98*H/2-minyloc;
      if smallDiam >= 0.25
        widthX = prams.Dpostx/2;
        centy = minyloc + smallDiam/2;
        sxwalls = widthX*cos(-tobs)+centx;
        sywalls = smallDiam/2*sin(-tobs)+centy;
        sidx = isempty(find(sywalls >= 0.99*H/2,1));
        if sidx
          XwallsInt(:,jdx) = [sxwalls;sywalls];
          jdx = jdx + 1;
        end
      end % smallDiam
    else
      maxyloc = max(ywall);
      smallDiam = maxyloc+0.98*H/2;
      if smallDiam >= 0.25
        widthX = prams.Dpostx/2;
        centy = maxyloc - smallDiam/2;
        sxwalls = widthX*cos(-tobs)+centx;
        sywalls = smallDiam/2*sin(-tobs)+centy;
        sidx = isempty(find(sywalls <= -0.99*H/2,1));
        if sidx
          XwallsInt(:,jdx) = [sxwalls;sywalls];
          jdx = jdx + 1;
        end % sidx
      end % smallDiam
    end % find(ywall<0)
  end %idx1&&idx2
end % iwall  
  
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
