clear; clc;
addpath ../src/
addpath ../examples/
addpath ./shannets/
addpath ./shannets/ves_fft_models/

pe = pyenv(Version="/opt/anaconda3/envs/mattorch/bin/python");

pathofDocument = fileparts(which('Net_ves_relax_midfat.py'));
if count(py.sys.path,pathofDocument) == 0
    insert(py.sys.path,int32(0),pathofDocument);
end

pathofDocument = fileparts(which('Net_ves_adv_fft.py'));
if count(py.sys.path,pathofDocument) == 0
    insert(py.sys.path,int32(0),pathofDocument);
end

pathofDocument = fileparts(which('ves_fft_mode2.pth'));
if count(py.sys.path,pathofDocument) == 0
    insert(py.sys.path,int32(0),pathofDocument);
end



% FLAGS
%-------------------------------------------------------------------------
prams.bgFlow = 'parabolic'; % 'shear','tayGreen','relax','parabolic'
prams.speed = 12000; % 500-3000 for shear, 70 for rotation, 100-400 for parabolic 
iplot = 0;
% PARAMETERS, TOOLS
%-------------------------------------------------------------------------
errTol = 1e-2;
maxDt = 1e-6; % dt = 1.28e-3,1e-3, 1.6e-4, 1e-5, 1e-6

if prams.speed == 100; prams.Th = 20; end; % Ca = 2.5
if prams.speed == 200; prams.Th = 10; end; % Ca = 5
if prams.speed == 400; prams.Th = 5; end; % Ca = 10
if prams.speed == 8000; prams.Th = 0.25; end; % Ca = 200
if prams.speed == 12000; prams.Th = 0.15; end; % Ca = 267
if prams.speed == 16000; prams.Th = 0.15; end;
if prams.speed == 24000; prams.Th = 0.05; end;
if prams.speed == 30000; prams.Th = 0.05; end;

% prams.Th = 0.05; % time horizon
prams.N = 128; % num. points for true solve in DNN scheme
prams.nv = 1; %(24 for VF = 0.1, 47 for VF = 0.2) num. of vesicles
prams.fmm = false; % use FMM for ves2ves
prams.fmmDLP = false; % use FMM for ves2walls
prams.kappa = 1;
prams.dt = maxDt; % time step size
prams.Nbd = 0;
prams.nvbd = 0;
prams.interpOrder = 1;
oc = curve;
Th = prams.Th; N = prams.N; nv = prams.nv; dt = prams.dt; 
bgFlow = prams.bgFlow; speed = prams.speed;

% net parameters
Nnet = 128; % num. points
nComp = 16; % number of components for networks except relaxation problem
nCompRelax = 32; % number of PCA components for relaxation problem
% # of modes for M*vinf's network, inv(DivGT)*Div*vinf uses the same
nVelModes = 24; activeModes = [(1:nVelModes/2)';(256-nVelModes/2+1:256)'];
% load PCA basis vectors and column mean (evects, colMeans) for Nnet = 256
% load necessaryMatFiles/pcaCoeffsBasis1step.mat
load necessaryMatFiles/pcaBasisNewest.mat


%-------------------------------------------------------------------------
disp(['Flow: ' prams.bgFlow ', N = ' num2str(N) ', nv = ' num2str(nv) ...
    ', Th = ' num2str(Th)])
%-------------------------------------------------------------------------

% VESICLES and WALLS:
% -------------------------------------------------------------------------
% initType = 2; % 1: VesID, 2: vesShape, 3: initialVesFile, 4: randomly fill
% vesID = 88201;  % 88201 from the library
% vesShape = 'ellipse'; %'ellipse'; % 'ellipse' or 'curly' 
% initialVesFile = []; %'VF35_81VesIC'; % from another run % probVF20IC: fails in 2it
% cent = [0; 0.065]; % centers of vesicles, if empty, randomly assigned 
% thet = []; % angular positions of vesicles, if empty, randomly assigned
% IA = [pi/2]; % initial inclination angles of vesicles
% irandInit = 0; % perturb center, angle and IA of vesicles?
% volFrac = 0;
% [X,area0,len0] = initializeVesiclesAndWalls(vesID,...
%     vesShape,initialVesFile,cent,thet,IA,initType,bgFlow,N,nv,oc,[]);


X0 = oc.initConfig(N,'ellipse');
[~,~,len] = oc.geomProp(X0);
X0 = X0./len;
IA = pi/2;
cent = [0; 0.065];
X = zeros(size(X0));
X(1:N) = cos(IA) * X0(1:N) - ...
      sin(IA) * X0(N+1:2*N) + cent(1);
X(N+1:2*N) = sin(IA) * X0(1:N) +  ...
      cos(IA) * X0(N+1:2*N) + cent(2);

[~,area0,len0] = oc.geomProp(X);

% [/2),X(end/2+1:end),vinf(1:end/2),vinf(end/2+1:end))



% -------------------------------------------------------------------------

solveType = 'DNN';
fileName = ['./output/poisDNNnewSingVes_speed' num2str(prams.speed) '_newNet_trueAdv_noSplit.bin'];
fid = fopen(fileName,'w');
output = [N;nv];
fwrite(fid,output,'double');
x = X(1:end/2,:); y = X(end/2+1:end,:);
output = [x(:); y(:)];
fwrite(fid,output,'double');
fclose(fid);

% BUILD DNN CLASS
% -------------------------------------------------------------------------
dnn = dnnToolsSingleVesPro(X,prams);

% LOAD NORMALIZATION PARAMETERS
load ./shannets/ves_fft_in_param.mat
load ./shannets/ves_fft_out_param.mat
dnn.torchAdvInNorm = in_param;
dnn.torchAdvOutNorm = out_param;

tt = dnn.tt; dnn.oc = oc; 
dnn.nComp = nComp;  dnn.nCompRelax = nCompRelax;
% LOAD FFT Based Network for M*Vinf 
dnn = loadFCnet4AdvectFiles(nVelModes,activeModes,dnn);
% save PCA matrices 
dnn.colMeans = colMeans; dnn.evects = evects; 
% -------------------------------------------------------------------------

% INITIALLY TAKE SMALL TIME STEPS WITH IMPLICIT SOLVER
% ------------------------------------------------------------------------
% necessary to find correct initial tension, density, rotlet and stokeslet
tt.dt = maxDt; sig = zeros(N,nv); eta = []; RS = [];
% for iter = 1 : 2
%   vesicle = capsules(X,[],[],prams.kappa,ones(nv,1),1); vesicle.setUpRate();
%   [X,sig,eta,RS] = tt.timeStepSimple(X,sig,[],[],ones(nv,1),[],vesicle);
% end
tt.dt = dt;
% ------------------------------------------------------------------------

% INITIALIZE MATRICES AND COUNTERS
% ------------------------------------------------------------------------
time = [0];
Xhist = X; sigStore = sig; 
errALPred = 0;
ncountCNN = 0;
ncountExct = 0;
% ------------------------------------------------------------------------
writeData(fileName,Xhist,sigStore,time(end),ncountCNN,ncountExct);

% TIME STEPPING
it = 1;
cx = []; cy = [];
while time(end) < prams.Th
  disp('********************************************') 
  disp([num2str(it) 'th time step, time: ' num2str(time(it))])
  
  
  
  
  disp('Taking a step with DNNs...');  tStart = tic;    
  Xnew = dnn.DNNsolveTorchNoSplit(Xhist);
  ncountCNN = ncountCNN + 1;
  
  [xIntersect,~,~] = oc.selfintersect(Xnew);
  if ~isempty(xIntersect); 
  disp('New vesicle shape is self-intersecting!!!'); 
  else
  % AREA-LENGTH CORRECTION
  disp('Correcting area-length errors...')
  [Xnew2,ifail] = oc.correctAreaAndLength2(Xnew,area0,len0);
  if ifail
    disp('Error in AL cannot be corrected!!!')
  else
    Xnew = oc.alignCenterAngle(Xnew,Xnew2);
  end
  end
  
  
  [~,area,len] = oc.geomProp(Xnew);
  errArea = max(abs(area-area0)./area0); errLen = max(abs(len-len0)./len0);
  
  it = it + 1;
  Xhist = Xnew;
  errALPred(it) = max(errArea,errLen);
  time(it) = time(it-1) + prams.dt;  
  

  cx = [cx; mean(Xhist(1:end/2))];
  cy = [cy; mean(Xhist(end/2+1:end))];

  disp(['Error in area and length: ' num2str(max(errArea, errLen))])   
  disp('********************************************') 
  disp(' ')
  
  if rem(it,10) == 0
    writeData(fileName,Xhist,sigStore,time(end),ncountCNN,ncountExct);  
    figure(2); clf;
    plot(cx, cy, 'linewidth',2)
    axis square
    grid
    xlabel('center in x')
    ylabel('center in y')
    title(['Time: ' num2str(time(it))])
  end

  disp(['took ' num2str(toc(tStart)) ' seconds.'])
 
  
  if iplot
  figure(1);clf;
  hold on;
  x = [Xhist(1:end/2); Xhist(1)];
  y = [Xhist(1+end/2:end); Xhist(end/2+1)];
  plot(x,y,'r','linewidth',2)
  plot(Xhist(1), Xhist(end/2+1),'o','markerfacecolor','r','markersize',8)
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
  load ./necessaryMatFiles/X100KinitShapes.mat  
  X0 = Xstore(:,vesID);
  X = [interpft(X0(1:end/2),N);interpft(X0(end/2+1:end),N)];
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
function writeData(filename,X,sigma,time,ncountNN,ncountExact)
x = X(1:end/2,:);
y = X(end/2+1:end,:);
output = [time;ncountNN;ncountExact;x(:);y(:);sigma(:)];

fid = fopen(filename,'a');
fwrite(fid,output,'double');
fclose(fid);


end
