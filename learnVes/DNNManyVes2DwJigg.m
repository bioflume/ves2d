clear; clc;
addpath ../src/
addpath ../examples/
disp('Multiple vesicles with background fluid, NNA scheme')

% FLAGS
%-------------------------------------------------------------------------
variableKbDt = ~true;
prams.bgFlow = 'couette'; % 'rotation' or 'couette' (w/ solid boundaries)
prams.speed = 100; % 70 for rotation, 100 for couette 
iDNNlike = 0; % whether solve exactly DNNlike or use DNNs
iplot = 0;
iJiggle = 1;
iuseNear = 0; % use wrong near-interactions (if =0, then neglect near-int)
iTrueNear = 1;
irepulsion = 1; % use repulsion, this is possible if iTrueNear == 1

% PARAMETERS, TOOLS
%-------------------------------------------------------------------------
exactSolveFreq = 100; % solve exactly at every [] time steps
prams.Th = 1.5; % time horizon
prams.N = 64; % num. points for true solve in DNN scheme
prams.nv = 81; %(24 for VF = 0.1, 47 for VF = 0.2) num. of vesicles
prams.fmm = ~false; % use FMM for ves2ves
prams.fmmDLP = ~false; % use FMM for ves2walls
prams.kappa = 1;

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
nCompRelax = 16; % number of PCA components for relaxation problem
% # of modes for M*vinf's network, inv(DivGT)*Div*vinf uses the same
nVelModes = 24; activeModes = [(1:nVelModes/2)';(Nnet-nVelModes/2+1:Nnet)'];
% number of modes for inv(DivGT)*(DivGB)x
nTenModes = 32; tenPredModes = [(1:nTenModes/2)';(Nnet-nTenModes/2+1:Nnet)'];
% load PCA basis vectors and column mean (evects, colMeans) for Nnet = 256
%load necessaryMatFiles/pcaCoeffsBasis1step.mat
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
initialVesFile = 'VF35_81VesIC'; % from another run % probVF20IC: fails in 2it
cent = []; % centers of vesicles, if empty, randomly assigned 
thet = []; % angular positions of vesicles, if empty, randomly assigned
IA = []; % initial inclination angles of vesicles
irandInit = 1; % perturb center, angle and IA of vesicles?
volFrac = 0.3;
[X,Xwalls,area0,len0] = initializeVesiclesAndWalls(vesID,...
    vesShape,initialVesFile,cent,thet,IA,initType,bgFlow,N,nv,Nbd,nvbd,irandInit,oc,[],[]);
% -------------------------------------------------------------------------

% BUILD DNN CLASS
% -------------------------------------------------------------------------
dnn = dnnToolsManyVesNoLoop(X,Xwalls,prams);
tt = dnn.tt; dnn.oc = oc; dnn.useNear = iuseNear; dnn.useTrueNear = iTrueNear;
dnn.repulsion = irepulsion; dnn.variableKbDt = variableKbDt;
if irepulsion
  vesicle = capsules(X(:,1),[],[],prams.kappa,1,true);
  vesicle.setUpRate();
  dnn.minDist = 0.4;
  dnn.repStrength = vesicle.repulStrengthScale(dnn.minDist,tt,prams.speed);
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
if ~variableKbDt
dnn = loadPCAnet4RelaxFiles(dnn);
else
dnn = loadAllPCAnets4Relax(dnn);
end
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
tt.dt = 1E-5; sig = zeros(N,nv); eta = zeros(2*Nbd,nvbd); RS = zeros(3,nvbd);
for iter = 1 : 2
  vesicle = capsules(X,[],[],prams.kappa,ones(nv,1),1); vesicle.setUpRate();
  [X,sig,eta,RS] = tt.timeStepSimple(X,sig,eta,RS,ones(nv,1),dnn.walls,vesicle);
end
tt.dt = dt;
% ------------------------------------------------------------------------

% INITIALIZE MATRICES AND COUNTERS
% ------------------------------------------------------------------------
time = (0:dt:Th)'; ntime = numel(time);
Xhist = zeros(2*N,nv,ntime); Xhist(:,:,1) = X; 
sigStore = zeros(N,nv,ntime); sigStore(:,:,1) = sig;
errALPred = zeros(ntime,1);
etaStore = zeros(2*Nbd,nvbd,ntime); RSstore = zeros(3,nvbd,ntime);
if Nbd ~= 0
  etaStore(:,:,1) = eta; RSstore(:,:,1) = RS;
end

ncountCNN = 0;
ncountExct = 0;
% ------------------------------------------------------------------------
if iDNNlike
  solveType = 'DNNlike';
else
  solveType = 'DNN';
end
if initType == 1
fileName = ['./output/nv' num2str(nv) 'N' num2str(N) ...
    solveType 'VesID' num2str(vesID) '_bgFlow' bgFlow '_speed' num2str(speed)];
elseif initType == 2
fileName = ['./output/OldNets_nv' num2str(nv) 'N' num2str(N) ...
    solveType 'VesShape' vesShape '_bgFlow' bgFlow '_speed' num2str(speed)];
elseif initType == 3
fileName = ['./output/nv' num2str(nv) 'N' num2str(N) ...
    solveType 'TrueNear_Kb1Dt1E4_VF35_bgFlow' bgFlow '_speed' num2str(speed)];
end

% TIME STEPPING
for it = 2 : ntime
  disp('********************************************') 
  disp([num2str(it) 'th (/' num2str(ntime) ...
    ') time step, time: ' num2str(time(it))])

  
  if rem(it,exactSolveFreq) == 0       
    % EXACT SOLVE  
    disp('Taking an exact time step...'); tStart = tic;  
    ncountExct = ncountExct + 1;
    vesicle = capsules(Xhist(:,:,it-1),[],[],prams.kappa,ones(nv,1),1);
    vesicle.setUpRate();
    [Xhist(:,:,it),sigStore(:,:,it),etaStore(:,:,it),RSstore(:,:,it),...
        iter,iflag] = tt.timeStepSimple(Xhist(:,:,it-1),...
        sigStore(:,:,it-1),etaStore(:,:,it-1),RSstore(:,:,it-1),ones(nv,1),...
        dnn.walls,vesicle);
    
  else 
    % NETWORK SOLVE 
    ncountCNN = ncountCNN+1;
    if iDNNlike
      disp('Taking a step with DNN-like exact scheme...');  tStart = tic;
      [Xhist(:,:,it),sigStore(:,:,it),etaStore(:,:,it),...
        RSstore(:,:,it)] = dnn.DNNlikeExactSolve(Xhist(:,:,it-1),...
        sigStore(:,:,it-1),etaStore(:,:,it-1),RSstore(:,:,it-1));
    else
      disp('Taking a step with DNNs...');  tStart = tic;
      [Xhist(:,:,it),sigStore(:,:,it),etaStore(:,:,it),...
          RSstore(:,:,it)] = dnn.DNNsolve(Xhist(:,:,it-1),...
          sigStore(:,:,it-1),etaStore(:,:,it-1),RSstore(:,:,it-1),Nnet);
    end % if iDNNlike

    % check if there is any self-intersecting shape
    for k = 1 : nv
      [xIntersect,~,~] = oc.selfintersect(Xhist(:,k,it));
      if ~isempty(xIntersect)
        disp('New vesicle shape is self-intersecting!!!')
      end
    end

    if ~isempty(xIntersect) 
      disp('Taking an exact time step...')
      ncountCNN = ncountCNN-1;
      ncountExct = ncountExct+1;
      vesicle = capsules(Xhist(:,:,it-1),[],[],prams.kappa,ones(nv,1),1);
      vesicle.setUpRate();
      [Xhist(:,:,it),sigStore(:,:,it),etaStore(:,:,it),RSstore(:,:,it),...
          iter,iflag] = tt.timeStepSimple(Xhist(:,:,it-1),...
          sigStore(:,:,it-1),etaStore(:,:,it-1),RSstore(:,:,it-1),ones(nv,1),...
          dnn.walls,vesicle);
    end % self-intersection occurs
    
  end % if rem(it,10) == 0
  
  if 0; XafterSolve(:,:,it) = Xhist(:,:,it); end;

  % JIGGLING
  if iJiggle % jiggle vesicles pointwise if there is a near-collision
    Xhist(:,:,it) = oc.fixCollisionsZinc(Xhist(:,:,it),Xwalls);
    
    if 0; XafterJigg(:,:,it) = Xhist(:,:,it); end;

    % Equally distribute points in arc-length b/c jiggling messes up points
    disp('Equally distributing points on arc-length...') 
    Xiter = Xhist(:,:,it);
    for iter = 1 : 5
      [Xiter,~,~] = oc.redistributeArcLength(Xiter);
    end
    
    
    if 0; XafterDist1(:,:,it) = Xiter; end;
    

    % Fix misalignment in center and angle due to reparametrization
    Xhist(:,:,it) = oc.alignCenterAngle(Xhist(:,:,it),Xiter);
    
    
    if 0; XafterAlignAfterDist1(:,:,it) = Xhist(:,:,it); end;

  end
  
  % AREA-LENGTH CORRECTION
  disp('Correcting area-length errors...')
  [Xnew,ifail] = oc.correctAreaAndLength2(Xhist(:,:,it),area0,len0);
  
  if 0; XafterCorr(:,:,it) = Xnew; end;

  if ifail
    disp('Error in AL cannot be corrected!!!')
  end
  Xhist(:,:,it) = oc.alignCenterAngle(Xhist(:,:,it),Xnew);

  if 0; XafterAlignAfterCorr(:,:,it) = Xhist(:,:,it); end; 

  % check if shape is intersecting
  for k = 1 : nv
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

  if 0; XafterDist2(:,:,it) = Xiter; end;
  
    % Fix misalignment in center and angle due to reparametrization
  Xhist(:,:,it) = oc.alignCenterAngle(Xhist(:,:,it),Xiter);
  
  if 0; XalignAfterDist2(:,:,it) = Xhist(:,:,it); end;
 

  if rem(it,100)
    save(fileName,'Xhist','dt','time','Xwalls','errALPred','it')
  end

  disp(['took ' num2str(toc(tStart)) ' seconds.'])
  
  % Compute error in area and length
  [~,area,len] = oc.geomProp(Xhist(:,:,it));
  errArea = max(abs(area-area0)./area0); errLen = max(abs(len-len0)./len0);
  errALPred(it) = max(errArea,errLen);
  disp(['Error in area and length: ' num2str(errALPred(it))])   
  disp('********************************************') 
  disp(' ')
  
  if iplot
  figure(1);clf;
  plot(Xwalls(1:end/2,:),Xwalls(end/2+1:end,:),'k')
  hold on;
  plot(Xhist(1:end/2,:,it),Xhist(end/2+1:end,:,it),'r','linewidth',2)
  axis equal
  pause(0.1)
  end
end

% Save data to a mat-file:
save(fileName,'Xhist','dt','time','sigStore','etaStore',...
'errALPred','vesID','ncountCNN','ncountExct','Xwalls')
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
function dnn = loadPCAnet4RelaxFiles(dnn)
% network for solving relaxation problem for coupled tension and X
    
if dnn.kappa == 1 % KAPPA = 1 NETWORKS (Dts = 1e-5, 5e-5, 1e-4, 5e-4, 1e-3)
if dnn.dt == 1e-5
  file = './networks/fcPCArelaxN256Dt1E5nModes';
elseif dnn.dt == 5e-5
  file = './networks/fcPCArelaxN256Dt5E5nModes';    
elseif dnn.dt == 1e-4
  file = './networks/fcPCArelaxN256Dt1E4Kb1nModes';
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
    
    
end % if dnn.kappa = 1

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dnn = loadAllPCAnets4Relax(dnn)

% list of networks for different Kb*Dt values
dnn.KbDts = [5e-7;1e-6;5e-6;1e-5;5e-5;1e-4;2.5E-4;5e-4;7.5E-4;1e-3;2.5E-3];

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
flowKbDt = dnn.kappa * dnn.dt;

[~, idx] = min(abs(dnn.KbDts-flowKbDt));
    
if idx == 1 || idx == numel(dnn.KbDts)
  disp('Extrapolation needed for the given bending stiffness and Dt, stop!')
  pause
else
  % Choose 5 networks and we do 5th order Lagrange Interpolation 
  if idx == numel(dnn.KbDts)-1
    whichNets = (idx-3:idx+1);
  elseif idx == 2
    whichNets = (idx-1:idx+3);
  else
    whichNets = (idx-2:idx+2);
  end
%   whichNets = (1:numel(dnn.KbDts))';
end

% Load the networks needed for Lagrange interpolation
for k = 1 : numel(whichNets)
  load([files{whichNets(k)} '1to16_fcXlarge_tstepFCNet_w1step.mat'])
  dnn.bendNets{k,1} = net; 
  dnn.muChan_bend(k,1) = muChan1; 
  dnn.sdevChan_bend(k,1) = sdevChan1; 
  dnn.scale_bend(k,1) = scale; 
  dnn.offset_bend(k,1) = offset;
  if dnn.nCompRelax > 16
  load([files{whichNets(k)} '17to32_fcXlarge_tstepFCNet_w1step.mat'])
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
