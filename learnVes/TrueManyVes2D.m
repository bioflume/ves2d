clear; clc;
addpath ../src/
addpath ../examples/
disp('Many vesicles with background fluid, Implicit scheme')

% FLAGS
%-------------------------------------------------------------------------
prams.bgFlow = 'couette'; % 'rotation' or 'couette' (confined)
prams.speed = 100; % 70 for rotation, 100 for couette
iplot = 0;
irandInit = 0;
iequalDist = 0;
ireparam = 1;
irepulsion = 0;
% PARAMETERS, TOOLS
%-------------------------------------------------------------------------
prams.Th = 1.5; % time horizon 
prams.N =  32; % num. points
prams.nv = 81; % number of vesicles
prams.fmm = true; % use FMM for ves2ves
prams.fmmDLP = true; % use FMM for ves2walls
prams.kappa = 5;

prams.dt = 1e-4; % time step size
if strcmp(prams.bgFlow,'couette')
  prams.Nbd = 256;
  prams.nvbd = 2;
else
  prams.Nbd = 0;
  prams.nvbd = 0;
end
oc = curve;
Th = prams.Th; N = prams.N; nv = prams.nv; dt = prams.dt; 
Nbd = prams.Nbd; nvbd = prams.nvbd; bgFlow = prams.bgFlow; speed = prams.speed;
disp(['Flow: ' prams.bgFlow ', N = ' num2str(N) ', nv = ' num2str(nv) ', Th = ' num2str(Th)])
%-------------------------------------------------------------------------

% VESICLES and WALLS:
% -------------------------------------------------------------------------
initType = 3; % 1: VesID, 2: vesShape, 3: initialVesFile
vesID = 88201;  % 88201 from the library
vesShape = 'ellipse'; % 'ellipse' or 'curly' 
initialVesFile = 'VF35_81VesIC'; % from another run
cent = []; % centers of vesicles, if empty, randomly assigned 
thet = []; % angular positions of vesicles, if empty, randomly assigned
IA = []; % initial inclination angles of vesicles

[X,Xwalls,area0,len0] = initializeVesiclesAndWalls(vesID,...
    vesShape,initialVesFile,cent,thet,IA,initType,bgFlow,N,nv,Nbd,nvbd,irandInit,oc);
% figure(1);clf;
% plot(Xwalls(1:end/2,:),Xwalls(end/2+1:end,:),'k')
% hold on
% plot(X(1:end/2,:),X(end/2+1:end,:),'r')
% axis equal
% pause;
% -------------------------------------------------------------------------

% BUILD DNN CLASS
% -------------------------------------------------------------------------
dnn = dnnToolsManyVesNoLoop(X,Xwalls,prams);
tt = dnn.tt;
tt.repulsion = irepulsion;
if irepulsion
  vesicle = capsules(X(:,1),[],[],prams.kappa,1,true);
  vesicle.setUpRate();
  tt.minDist = 0.4;
  tt.repStrength = vesicle.repulStrengthScale(tt.minDist,tt,prams.speed);
end
% -------------------------------------------------------------------------

% INITIALIZE MATRICES AND COUNTERS
% ------------------------------------------------------------------------
timeTrue = (0:dt:Th)'; ntime = numel(timeTrue);
XhistTrue = zeros(2*N,nv,ntime);
XhistTrue(:,:,1) = X;
sigStore = zeros(N,nv,ntime);
errALTrue = zeros(ntime);
etaStore = zeros(2*Nbd,nvbd,ntime);
RSstore = zeros(3,nvbd,ntime);
% ------------------------------------------------------------------------
if initType == 1
fileName = ['./output/nv' num2str(nv) 'N' num2str(N) 'TrueVesID' num2str(vesID) '_bgFlow' bgFlow '_speed' num2str(speed)];
elseif initType == 2
fileName =  ['./output/nv' num2str(nv) 'N' num2str(N) 'TrueVesShape' vesShape '_bgFlow' bgFlow '_speed' num2str(speed)];
elseif initType == 3
fileName =  ['./output/TrueN32_VF35Kb5Dt1e4_bgFlow' bgFlow '_speed' num2str(speed)];    
end

% TIME STEPPING
for it = 2 : ntime
  disp('********************************************') 
  disp([num2str(it) 'th (/' num2str(ntime) ...
    ') time step, time: ' num2str(timeTrue(it))])
  
  tStart = tic;
  disp('---------------------------')
 
  % SOLVE USING IMPLICIT TIME STEPPING
  disp('Solving implicitly (the most recent scheme)...')  
  vesicle = capsules(XhistTrue(:,:,it-1),[],[],prams.kappa,ones(nv,1),1);
  vesicle.setUpRate();
  [XhistTrue(:,:,it),sigStore(:,:,it),etaStore(:,:,it),RSstore(:,:,it),...
      iter,iflag] = tt.timeStepSimple(XhistTrue(:,:,it-1),...
      sigStore(:,:,it-1),etaStore(:,:,it-1),RSstore(:,:,it-1),ones(nv,1),...
      dnn.walls,vesicle);

  % AREA-LENGTH CORRECTION
  [Xnew,ifail] = oc.correctAreaAndLength2(XhistTrue(:,:,it),area0,len0);
  if ifail
    disp('Error in AL cannot be corrected!!!')
  end
  XhistTrue(:,:,it) = oc.alignCenterAngle(XhistTrue(:,:,it),Xnew);

  % check if shape is intersecting
  for k = 1 : nv
    [xIntersect,~,~] = oc.selfintersect(XhistTrue(:,k,it));
    if ~isempty(xIntersect)
      disp('New vesicle shape is self-intersecting!!!')
    end
  end

  % Equally distribute points in arc-length or reparameterize
  if iequalDist 
  Xiter = XhistTrue(:,:,it);
  for iter = 1 : 5
    [Xiter,~,~] = oc.redistributeArcLength(Xiter);
  end
  % Fix misalignment in center and angle due to reparametrization
  XhistTrue(:,:,it) = oc.alignCenterAngle(XhistTrue(:,:,it),Xiter);
  elseif ireparam
  [Xiter,~] = oc.reparametrize(XhistTrue(:,:,it),XhistTrue(:,:,it)-XhistTrue(:,:,it-1),6,20);
  % Fix misalignment in center and angle due to reparametrization
  XhistTrue(:,:,it) = oc.alignCenterAngle(XhistTrue(:,:,it),Xiter);      
  end
  
  if rem(it,100) == 0
    save(fileName,'XhistTrue','Xwalls','dt','timeTrue','errALTrue','it')
  end

  disp(['took ' num2str(toc(tStart)) ' seconds.'])
  disp('---------------------------')    
  
  % Compute error in area and length
  [~,area,len] = oc.geomProp(XhistTrue(:,:,it));
  errArea = max(abs(area-area0)./area0); errLen = max(abs(len-len0)./len0);
  errALTrue(it) = max(errArea,errLen);
  disp(['Error in area and length: ' num2str(errALTrue(it))])   
  disp('********************************************') 
  disp(' ')
  
  if iplot
  figure(1);clf;
  plot(Xwalls(1:end/2,:),Xwalls(end/2+1:end,:),'k')
  hold on;
  plot(XhistTrue(1:end/2,:,it),XhistTrue(end/2+1:end,:,it),'r','linewidth',2)
  axis equal
  pause(0.1)
  end
end

% SAVE DATA 
save(fileName,'XhistTrue','Xwalls','dt','timeTrue','errALTrue','vesID')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Xwalls,area0,len0] = initializeVesiclesAndWalls(vesID,...
    vesShape,initialVesFile,cent,thet,IA,initType,bgFlow,N,nv,Nbd,nvbd,irandInit,oc)
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
end

if isempty(IA)
  IA = 2*pi*rand(nv,1);
end

if initType < 3
  [~,area0,len0] = oc.geomProp(X0);
  vesRadius = sqrt(area0(1)/pi); % 0.1291
  area0 = ones(1,nv)*area0; len0 = ones(1,nv)*len0;

  X = zeros(2*N,nv);
  % Move center for rotation flow
  if strcmp(bgFlow,'rotation')
    if isempty(thet)
      thet = 0:2*pi/nv:2*pi*(nv-1)/nv; 
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
      thet = 0:2*pi/nv:2*pi*(nv-1)/nv; 
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
end % if initType

Xwalls = zeros(2*Nbd,nvbd);
if strcmp(bgFlow,'couette')
  thet = (0:Nbd-1)'*2*pi/Nbd;
  Xwalls = [ [2.2*cos(thet); 2.2*sin(thet)] [cos(-thet); sin(-thet)] ];
end

end % initializeVesiclesAndWalls

