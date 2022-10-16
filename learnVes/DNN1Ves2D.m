clear; clc;
addpath ../src/
disp('Single vesicle with background fluid, resNet+Exact Solve')

% FLAGS
%-------------------------------------------------------------------------
bgFlow = 'relax'; % 'shear','tayGreen','relax','parabolic'
speed = 0 % 500-3000 for shear, 70 for rotation, 100-400 for parabolic 
iSolveTruth = 0; % also solve exactly
iAdvectWNN = 1; % use NN for advection
iRelaxWNN = 1; % use NN for relaxation
iFFTin2Advect = 0; % input to advection NN is FFT coeffs or PCA
iResNetRelax = 0; % use resNet for relaxation?

% PARAMETERS, TOOLS
%-------------------------------------------------------------------------
exactSolveFreq = 10; % solve exactly at every [] time steps
Th = 0.05; % time horizon
N = 256; % num. points
Nsolve = 16; % num. points for true solve in DNN scheme
dt = 1e-4; % time step size
kappa = 1;
nComp = 16; % number of components for relaxation problem's PCA network
nVelModes = 12; % number of modes for M*vinf's FFT based network
shapeModes = [2;256;4;255;6;3;5;254]; % input shape's modes
oc = curve;
op = poten(Nsolve);
dnn = dnnTools;
% load PCA basis vectors and column mean (evects, colMeans)
if N == 256
load necessaryMatFiles/pcaCoeffsBasis1step.mat
elseif N == 96
load necessaryMatFiles/pcaBasisN96for100K.mat
evects = evectsN96; colMeans = colMeansN96;
end
%-------------------------------------------------------------------------

% VESICLE:
% -------------------------------------------------------------------------
load ./necessaryMatFiles/errInPred100K_FCnetPCA.mat
errAllPCA = errPCA;
load /workspace/gokberk/X100KinitShapes.mat
vesID = [] 

if ~isempty(vesID)
  X = Xstore(:,vesID);
  X = [interpft(X(1:end/2),N);interpft(X(end/2+1:end),N)];
  %errFFT = ceil(errAllFFT(vesID)*100);
  errPCA = ceil(errAllPCA(vesID)*100);
  disp(['VesID: ' num2str(vesID) ', Error in prediction: (PCA) ' ...
      num2str(errPCA) '%'])
else
  initShape = 'curly'; %'curly', 'star','openStar'
  if strcmp(initShape,'ellipse');
    X = oc.ellipse(N,0.65);
  else
    X = oc.initConfig(N,initShape);
  end
  [~,area,len] = oc.geomProp(X);
  X = X/len;
end
[~,area0,len0] = oc.geomProp(X);
vesRadius = sqrt(area0/pi);
% Move center for Taylor-Green flow
if strcmp(bgFlow,'tayGreen')
X = [X(1:end/2)+pi/4; X(end/2+1:end)+pi/4];
end
% Move center for parabolic flow
if strcmp(bgFlow,'parabolic')
X = [X(1:end/2);X(end/2+1:end)+0.06];
end
% Move center for rotation flow
if strcmp(bgFlow,'rotation')
X = [X(1:end/2)+10*vesRadius;X(end/2+1:end)];
end
% -------------------------------------------------------------------------

% LOAD NETWORKS (SEE THE END OF THE CODE FOR THESE FUNCTIONS)
% -------------------------------------------------------------------------
if ~iResNetRelax
  % LOAD PCA Network for Relaxation Problem
  [PCAnets,muChan1_pca,sdevChan1_pca,scale_pca,offset_pca] = loadPCAnet4RelaxFiles(N,dt,kappa);
else
  % LOAD ResNet for Relaxation Problem  
  [PCAnets,muChan1_pca,sdevChan1_pca,scale_pca,offset_pca] = loadPCAResNet4RelaxFiles(nComp);  
end

% LOAD FFT Based Network for M*Vinf (input is FFT or PCA coeffs)
activeModes = [(1:nVelModes/2)';(N-nVelModes/2+1:N)'];
[FCnets,muChan1_fc,sdevChan1_fc,scale_fc,offset_fc,outputSize] = ...
    loadFCnet4AdvectFiles(iFFTin2Advect,nVelModes,activeModes,N);
% -------------------------------------------------------------------------

% BACKGROUND VELOCITY
%-------------------------------------------------------------------------
vinf = dnn.setBgFlow(bgFlow,speed);

% INITIALIZE MATRICES AND COUNTERS
% ------------------------------------------------------------------------
time = (0:dt:Th)';
Xhist = zeros(2*N,numel(time));
Xhist(:,1) = X;
errALPred = zeros(numel(time),1);
ncountCNN = 0;
ncountExct = 0;

XhistTrue = []; errALTrue = [];
if iSolveTruth
  XhistTrue = Xhist;
  errALTrue = errALPred;
end
% ------------------------------------------------------------------------

% TIME STEPPING
if ~isempty(vesID)
fileName = ['./output/N' num2str(Nsolve) 'xLDNNsimVesID' num2str(vesID) '_bgFlow' bgFlow '_speed' num2str(speed)];
else
fileName = ['./output/N' num2str(Nsolve) 'xLDNNsimVesShape' initShape '_bgFlow' bgFlow '_speed' num2str(speed)];
end
for it = 2 : numel(time)
  disp('********************************************') 
  disp([num2str(it) 'th (/' num2str(numel(time)) ...
    ') time step, time: ' num2str(time(it))])
  
  % SOLVE EXACTLY AND STORE SOLUTION TO COMPARE WITH NN+EXACT
  if iSolveTruth
    disp('---------------------------')
    disp('Solving for ground truth...') 
    tStart = tic;
    if N == 256
      Xin = [interpft(XhistTrue(1:end/2,it-1),96);...
        interpft(XhistTrue(end/2+1:end,it-1),96)];
      [Xout,errALTrue(it)] = dnn.exactlySolve(Xin,vinf,dt,area0,len0,oc,op);
      XhistTrue(:,it) = [interpft(Xout(1:end/2),256);...
        interpft(Xout(end/2+1:end),256)];
    else
      XhistTrue(:,it) = dnn.exactlySolve(XhistTrue(:,it-1),vinf,dt,area0,...
        len0,oc,op,kappa);
    end    
    disp(['took ' num2str(toc(tStart)) ' seconds.'])
    disp('---------------------------')
  end

  % NETWORK SOLUTION
  if rem(it,exactSolveFreq) == 0       
    disp('Taking an exact time step...'); tStart = tic;  
    ncountExct = ncountExct + 1;
    if N == 256
      Xin = [interpft(Xhist(1:end/2,it-1),Nsolve);...
          interpft(Xhist(end/2+1:end,it-1),Nsolve)]; 
      [Xout,errALPred(it)] = dnn.exactlySolve(Xin,vinf,dt,area0,len0,oc,op,kappa);
      Xhist(:,it) = [interpft(Xout(1:end/2),256);...
          interpft(Xout(end/2+1:end),256)];
    elseif N == 96
      [Xhist(:,it),errALPred(it)] = dnn.exactlySolve(Xhist(:,it-1),vinf,...
          dt,area0,len0,oc,op,kappa);
    end
  else
    disp('Taking a step with prediction...'); tStart = tic;
    ncountCNN = ncountCNN+1;
    
    % OPERATOR SPLITTING
    if iAdvectWNN
      % 1) TRANSLATION w/ NETWORK (USE PCA or FFT coeffs as INPUTS)  
      % Standardize the shape
      [X,scaling,rotate,~,sortIdx] = dnn.standardizationStep(Xhist(:,it-1),oc);  
      % Prepare input
      Xinput = dnn.prepareInputForNet(~iFFTin2Advect,X,colMeans,evects,...
        nComp,scale_fc,muChan1_fc,sdevChan1_fc,offset_fc,shapeModes);
      % Take a step
      Xnew = dnn.translateVinfwNN(Xinput,FCnets,activeModes,...
        vinf(Xhist(:,it-1)),dt,Xhist(:,it-1),outputSize,rotate,sortIdx);    
    else
      % 1) TRANSLATION w/o NETWORK
      Xnew = dnn.translateVinf(vinf(Xhist(:,it-1)),dt,Xhist(:,it-1),op);    
    end % if iAdvectWNN
   
    
    if iRelaxWNN == 1
      % 2) RELAXATION w/ NETWORK  
      % Standardize the shape
      [X,scaling,rotate,trans,sortIdx] = dnn.standardizationStep(Xnew,oc);
      % Prepare input
      Xinput = dnn.prepareInputForNet(~iFFTin2Advect,X,colMeans,evects,...
         nComp,scale_pca,muChan1_pca,sdevChan1_pca,offset_pca,shapeModes);
      % Take a step
      Xhist(:,it) = dnn.relaxWNN(Xinput,PCAnets,nComp,evects,colMeans,...
         scaling,rotate,trans,sortIdx,muChan1_pca,sdevChan1_pca,...
         offset_pca,scale_pca,iResNetRelax);
    else
      % 2) RELAXATION w/o NETWORK
      vesicle = capsules(Xnew,[],[],kappa,1,1); vesicle.setUpRate();
      Xhist(:,it) = dnn.relaxExactSolve(vesicle,zeros(size(Xnew)),dt,Xnew,op);
    end
    

    % AREA-LENGTH CORRECTION
    [Xnew,ifail] = oc.correctAreaAndLength2(Xhist(:,it),area0,len0);
    if ifail
      disp('Error in AL cannot be corrected!!!')
    end
    Xhist(:,it) = oc.alignCenterAngle(Xhist(:,it),Xnew);

    % check if shape is intersecting
    [xIntersect,~,~] = oc.selfintersect(Xhist(:,it));
    if ~isempty(xIntersect)
      disp('New vesicle shape is self-intersecting!!!')
    end

    if (ifail) || ~isempty(xIntersect) 
      disp('Taking an exact time step...')
      ncountCNN = ncountCNN-1;
      ncountExct = ncountExct+1;
      
      [Xhist(:,it),errALPred(it)] = dnn.exactlySolve(Xhist(:,it-1),vinf,dt,...
         area0,len0,oc,op,kappa);
    end % if area length correction fails or self-intersection occurs
  end % if rem(it,5) == 0
  
  % Equally distribute points in arc-length
  Xiter = Xhist(:,it);
  for iter = 1 : 5
    [Xiter,~,~] = oc.redistributeArcLength(Xiter);
  end
  % Fix misalignment in center and angle due to reparametrization
  Xhist(:,it) = oc.alignCenterAngle(Xhist(:,it),Xiter);

  disp(['took ' num2str(toc(tStart)) ' seconds.'])
  
  % Compute error in area and length
  [~,area,len] = oc.geomProp(Xhist(:,it));
  errALPred(it) = max(abs(area-area0)/area0,abs(len-len0)/len0);
  disp(['Error in area and length: ' num2str(errALPred(it))])   
  disp('********************************************') 
  disp(' ')
  if strcmp(bgFlow,'rotation')
    cx = mean(Xhist(1:end/2,it)); cy = mean(Xhist(end/2+1:end,it));
    if sqrt(cx^2+cy^2)<=3*vesRadius
      save(fileName,'Xhist','dt','time','errALPred','vesID','ncountCNN','ncountExct',...
      'errALTrue','XhistTrue','it')
      disp('Can kill the simulation, vesicle approaches to origin')
    end
  end
end


save(fileName,'Xhist','dt','time','XhistTrue',...
'errALPred','errALTrue','vesID','ncountCNN','ncountExct')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PCAnets,muChan1,sdevChan1,scale,offset] = loadPCAnet4RelaxFiles(N,dt,kappa)
% FULLY CONNECTED LAYERS, OUTPUT IS ALL COEFFs. 
muChan1 = []; sdevChan1 = []; scale = []; offset = [];

% normalized output
if kappa == 1
  if N == 256
    load('./networks/fcPCArelax_fcXlarge1to16_tstepFCNet_w1step.mat')
  elseif N == 96
    load('./networks/fcPCArelaxN96_fcXlarge_tstepFCNet_w1step.mat')
  end
elseif kappa == 1e-1 % ONLY N = 256
  load('./networks/fcPCArelaxN25Dt1E4Kb1E1nModes1to16_fcXlarge_tstepFCNet_w1step.mat') 
end

PCAnets = net;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FCnets,muChan1,sdevChan1,scale,offset,outputSize] = ...
  loadFCnet4AdvectFiles(iFFTin2Advect,nmodes,activeModes,N)

% initialize mu,sdev. for PCA net, we do not output any of them.
muChan1 = []; sdevChan1 = []; scale = []; offset = []; outputSize = [];

if N == 256
  if iFFTin2Advect 
    netFilesFolder = './networks/fcNetVelocityPredFilesFFTin/velPredFFTin_mode';
  else
    netFilesFolder = './networks/fcNetVelocityPredFilesPCAin/velPredPCAin_mode';
  end
else
  netFilesFolder = './networks/n96MtimesFFTNets/velPredPCAin_mode';    
end

for imode = 1 : nmodes
  pmode = activeModes(imode);
  load([netFilesFolder num2str(pmode) '_net_w1step.mat'])
  FCnets{imode} = net;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PCAnets,muChan1,sdevChan1,scale,offset] = ...
    loadPCAResNet4RelaxFiles(nComp)

muChan1 = []; sdevChan1 = []; scale = []; offset = [];

netFilesFolder = './networks/resNetTstepPCA1stepFiles/principComp';
fileEnd = '_w1step.mat';

for icoeff = 1 : nComp
  load([netFilesFolder num2str(icoeff) '_tstepNet' fileEnd])
  PCAnets{icoeff} = net;
end

end
