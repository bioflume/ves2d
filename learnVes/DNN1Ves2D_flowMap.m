clear; clc;
addpath ../src/
addpath ../examples/

% FLAGS
%-------------------------------------------------------------------------
variableKbDt = true;
interpOrder = 5; % must be an odd number
prams.bgFlow = 'parabolic'; % 'shear','tayGreen','relax','parabolic'
prams.speed = 200; % 500-3000 for shear, 70 for rotation, 100-400 for parabolic 
iDNNlike = 0; % whether solve exactly DNNlike or use DNNs
iplot = 0;
iJiggle = 0;
iuseNear = 0; % use wrong near-interactions (if =0, then neglect near-int)
iTrueNear = 0;
irepulsion = 0; % use repulsion, this is possible if iTrueNear == 1
toDown = 16;
% PARAMETERS, TOOLS
%-------------------------------------------------------------------------
exactSolveFreq = 0; % solve exactly at every [] time steps
errTol = 1e-2;
maxDt = 1e-4;
prams.Th = 2.5; % time horizon
prams.N = 256; % num. points for true solve in DNN scheme
prams.nv = 1; %(24 for VF = 0.1, 47 for VF = 0.2) num. of vesicles
prams.fmm = false; % use FMM for ves2ves
prams.fmmDLP = false; % use FMM for ves2walls
prams.kappa = 1;
prams.interpOrder = interpOrder;
prams.dt = 1E-4; % time step size
prams.Nbd = 0;
prams.nvbd = 0;

oc = curve;
Th = prams.Th; N = prams.N; nv = prams.nv; dt = prams.dt; 
bgFlow = prams.bgFlow; speed = prams.speed;

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
initType = 2; % 1: VesID, 2: vesShape, 3: initialVesFile, 4: randomly fill
vesID = 88201;  % 88201 from the library
vesShape = 'ellipse'; % 'ellipse' or 'curly' 
initialVesFile = []; %'VF35_81VesIC'; % from another run % probVF20IC: fails in 2it
cent = [0; 0.4]; % centers of vesicles, if empty, randomly assigned 
thet = []; % angular positions of vesicles, if empty, randomly assigned
IA = [-0.2]; % initial inclination angles of vesicles
irandInit = 0; % perturb center, angle and IA of vesicles?
volFrac = 0;
[X,area0,len0] = initializeVesiclesAndWalls(vesID,...
    vesShape,initialVesFile,cent,thet,IA,initType,bgFlow,N,nv,oc,[]);

% -------------------------------------------------------------------------

solveType = 'DNN';
fileName = ['./output/poisDNNnewSingVesInterp5_Vel200.bin'];
fid = fopen(fileName,'w');
output = [N;nv];
fwrite(fid,output,'double');
x = X(1:end/2,:); y = X(end/2+1:end,:);
output = [x(:); y(:)];
fwrite(fid,output,'double');
fclose(fid);

% BUILD DNN CLASS
% -------------------------------------------------------------------------
dnn = dnnToolsSingleVes(X,prams);
tt = dnn.tt; dnn.oc = oc; dnn.variableKbDt = variableKbDt;
if irepulsion
  vesicle = capsules(X(:,1),[],[],prams.kappa,1,true);
  vesicle.setUpRate();
  dnn.minDist = 0.4;
  dnn.repStrength = vesicle.repulStrengthScale(dnn.minDist,tt,prams.speed);
end    
% -------------------------------------------------------------------------

% LOAD NETWORKS (SEE THE END OF THE CODE FOR THESE FUNCTIONS)
% -------------------------------------------------------------------------
if ~iDNNlike
dnn.nComp = nComp;  dnn.nCompRelax = nCompRelax;
% LOAD PCA Network for Relaxation Problem
dnn = loadAllPCAnets4Relax(dnn);
% LOAD FFT Based Network for M*Vinf 
dnn = loadFCnet4AdvectFiles(nVelModes,activeModes,dnn);
% save PCA matrices 
dnn.colMeans = colMeans; dnn.evects = evects; 
end
% -------------------------------------------------------------------------

% INITIALLY TAKE SMALL TIME STEPS WITH IMPLICIT SOLVER
% ------------------------------------------------------------------------
% necessary to find correct initial tension, density, rotlet and stokeslet
tt.dt = maxDt; sig = zeros(N,nv); eta = []; RS = [];
for iter = 1 : 2
  vesicle = capsules(X,[],[],prams.kappa,ones(nv,1),1); vesicle.setUpRate();
  [X,sig,eta,RS] = tt.timeStepSimple(X,sig,[],[],ones(nv,1),[],vesicle);
end
tt.dt = dt;
% ------------------------------------------------------------------------

% INITIALIZE MATRICES AND COUNTERS
% ------------------------------------------------------------------------
time = [0];
Xhist = X; sigStore = sig; 
errALPred = 0;
sigGuess = zeros(N,nv);
ncountCNN = 0;
ncountExct = 0;
% ------------------------------------------------------------------------
writeData(fileName,Xhist,sigStore,time(end),ncountCNN,ncountExct);

% TIME STEPPING
it = 1;
while time(end) < prams.Th
  disp('********************************************') 
  disp([num2str(it) 'th time step, time: ' num2str(time(it))])
  
  shouldSolveExact = 0;
  
  if rem(it,exactSolveFreq) == 0
  disp('Taking an exact time step ...'); tStart = tic;
  vesicle = capsules(Xhist(:,it),sigGuess,[],prams.kappa,ones(1,1),1);
  vesicle.setUpRate();
  [Xnew,~,~,~,iter,iflag] = tt.timeStepSimple(Xhist,sigGuess,[],[],ones(1,1),[],vesicle);
  ncountExct = ncountExct + 1;
  else
  disp('Taking a step with DNNs...');  tStart = tic;    
  Xnew = dnn.DNNsolve(Xhist,Nnet);
  ncountCNN = ncountCNN + 1;
  end
  
  
  [xIntersect,~,~] = oc.selfintersect(Xnew);
  if ~isempty(xIntersect); 
  disp('New vesicle shape is self-intersecting!!!'); 
  shouldSolveExact = 1;
  else
  % AREA-LENGTH CORRECTION
  disp('Correcting area-length errors...')
  [Xnew2,ifail] = oc.correctAreaAndLength3(Xnew,area0,len0,toDown);
  if ifail
    shouldSolveExact = 1;  
    disp('Error in AL cannot be corrected!!!')
  else
    Xnew = oc.alignCenterAngle(Xnew,Xnew2);
  end
  end
  
  if shouldSolveExact
  disp('WE NEED TO SOLVE EXACT DUE TO NOT BEING ABLE CORRECT THE ERRORS')
  ncountCNN = ncountCNN - 1;
  ncountExct = ncountExct + 1; 
  vesicle = capsules(Xhist,sigGuess,[],prams.kappa,ones(1,1),1);
  vesicle.setUpRate();
  [Xnew,~,~,~,iter,iflag] = tt.timeStepSimple(Xhist,sigGuess,[],[],ones(1,1),[],vesicle);    
  [Xnew2,ifail] = oc.correctAreaAndLength3(Xnew,area0,len0,toDown);
  if ifail
    disp('Error in AL cannot be corrected EVEN AFTER EXACT SOLVE!!!')
    break;
  else
    Xnew = oc.alignCenterAngle(Xnew,Xnew2);
  end
  end
  
  [~,area,len] = oc.geomProp(Xnew);
  errArea = max(abs(area-area0)./area0); errLen = max(abs(len-len0)./len0);
  dtnew = prams.dt;
  if ~isempty(xIntersect)
    dtnew = 0.5 * prams.dt;
    if dtnew > maxDt; dtnew = maxDt; end;
  else
  if max(errArea, errLen) >= errTol
    dtnew = 0.5 * prams.dt; % Reject
    if dtnew > maxDt; dtnew = maxDt; end;
  elseif max(errArea, errLen) < 0.9 * errTol
    dtnew = 1.2 * prams.dt; % Accept
    if dtnew > maxDt; dtnew = maxDt; end;
    it = it + 1;
    Xhist = Xnew;
    errALPred(it) = max(errArea,errLen);
    time(it) = time(it-1) + dtnew;  
  else % still accept but do not change the time step size
    it = it + 1;
    Xhist = Xnew;
    errALPred(it) = max(errArea,errLen);
    time(it) = time(it-1) + dtnew;  
  end
  end
  disp(['Error in area and length: ' num2str(max(errArea, errLen))])   
  disp(['New time step size: ' num2str(dtnew)])
  disp('********************************************') 
  disp(' ')
  
  if rem(it,10) == 0
    writeData(fileName,Xhist,sigStore,time(end),ncountCNN,ncountExct);  
  end

  disp(['took ' num2str(toc(tStart)) ' seconds.'])
  
  if dtnew ~= prams.dt
  tt.dt = dtnew;
  prams.dt = dtnew;
  dnn.dt = dtnew;
  dnn.tt = tt;
  dnn = loadAllPCAnets4Relax(dnn);
  end
  
  if iplot
  figure(1);clf;
  hold on;
  x = [Xhist(1:end/2,:,it); Xhist(1,:,it)];
  y = [Xhist(1+end/2:end,:,it); Xhist(end/2+1,:,it)];
  plot(x,y,'r','linewidth',2)
  plot(Xhist(1,:,it), Xhist(end/2+1,:,it),'o','markerfacecolor','r','markersize',8)
  xlim([-1 1])
  ylim([-1 1])
  axis equal
  pause(0.1)
  end
end

% Save data to a mat-file:
writeData(fileName,Xhist,sigStore,time(end),ncountCNN,ncountExct);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,area0,len0] = initializeVesiclesAndWalls(vesID,...
    vesShape,initialVesFile,cent,thet,IA,initType,bgFlow,N,nv,oc,tt)

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
  X = zeros(size(X0));
  X(1:N) = cos(IA) * X0(1:N) - ...
      sin(IA) * X0(N+1:2*N) + cent(1);
  X(N+1:2*N) = sin(IA) * X0(1:N) +  ...
      cos(IA) * X0(N+1:2*N) + cent(2);
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
  X0(1:end/2) = X0(1:end/2) + cent(1);
  X0(end/2+1:end) = X0(end/2+1:end) + cent(2);
end

if isempty(IA)
  if irandInit
    IA = 2*pi*rand(nv,1);
  else
    IA = zeros(nv,1);
  end
end

[ra,area0,len0] = oc.geomProp(X);



end % initializeVesiclesAndWalls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dnn = loadAllPCAnets4Relax(dnn)
interpOrder = dnn.interpOrder;
% list of networks for different Kb*Dt values
dnn.KbDts = 2.^(0:10)' * 1E-6;
files = [];
fileName = './networks2/fcPCArelaxN256Dt';
for ifl = 1 : numel(dnn.KbDts)
  files{ifl} = [fileName num2str(dnn.KbDts(ifl)) 'Kb1nModes'];
end

% current flow's Kb*Dt
flowKbDt = dnn.kappa * dnn.dt;

[~, idx] = min(abs(dnn.KbDts-flowKbDt));
if idx == 1 || idx == numel(dnn.KbDts)
  disp('Extrapolation needed for the given bending stiffness and Dt, stop!')
  pause
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
  if 1%dnn.KbDts(k) >= 5E-4
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
  if 1%dnn.KbDts(k) >= 5E-4
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
function writeData(filename,X,sigma,time,ncountNN,ncountExact)
x = X(1:end/2,:);
y = X(end/2+1:end,:);
output = [time;ncountNN;ncountExact;x(:);y(:);sigma(:)];

fid = fopen(filename,'a');
fwrite(fid,output,'double');
fclose(fid);


end
