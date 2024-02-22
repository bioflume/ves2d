clear; clc;
addpath ../src/
disp('Single vesicle with background fluid, Exact Solve')

% FLAGS
%-------------------------------------------------------------------------
bgFlow = 'parabolic'; % 'shear','tayGreen','relax','parabolic','rotation'
speed = 8000 % 500-3000 for shear, 70 for rotation, 100-400 for parabolic 
iSplit = 0; % whether split operators or not
iequalDist = 0;
ireparam = 1;
kappa = 1;
% PARAMETERS, TOOLS
%-------------------------------------------------------------------------
Th = 10; % time horizon 
N =  128 % num. points
dt = 5e-4; % time step size
oc = curve;
op = poten(N);
dnn = dnnTools;
%-------------------------------------------------------------------------

% VESICLE:
% -------------------------------------------------------------------------
% load ./necessaryMatFiles/errInPred100K_FCnetPCA.mat
% errAllPCA = errPCA;
errAllPCA = [];
% load /workspace/gokberk/X100KinitShapes.mat
vesID = []; %88201;

if ~isempty(vesID)
  X = Xstore(:,vesID);
  X = [interpft(X(1:end/2),N);interpft(X(end/2+1:end),N)];
  %errFFT = ceil(errAllFFT(vesID)*100);
  errPCA = ceil(errAllPCA(vesID)*100);
  disp(['VesID: ' num2str(vesID) ', Error in prediction: (PCA) ' ...
      num2str(errPCA) '% '])
else
  initShape = 'ellipse';
  X = oc.initConfig(N,initShape);
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
% X = [X(1:end/2);X(end/2+1:end)+0.06];
X0 = X;
IA = pi/2;
cent = [0; 0.065];
X(1:N) = cos(IA) * X0(1:N) - ...
      sin(IA) * X0(N+1:2*N) + cent(1);
X(N+1:2*N) = sin(IA) * X0(1:N) +  ...
      cos(IA) * X0(N+1:2*N) + cent(2);
end
% Move center for rotation flow
if strcmp(bgFlow,'rotation')
X = [X(1:end/2)+10*vesRadius;X(end/2+1:end)];
end
% -------------------------------------------------------------------------

% BACKGROUND VELOCITY
%-------------------------------------------------------------------------
vinf = dnn.setBgFlow(bgFlow,speed);

% INITIALIZE MATRICES AND COUNTERS
% ------------------------------------------------------------------------
timeTrue = (0:dt:Th)';
XhistTrue = zeros(2*N,numel(timeTrue));
XhistTrue(:,1) = X;
errALTrue = zeros(numel(timeTrue),1);

% ------------------------------------------------------------------------
if ~isempty(vesID)
fileName = ['./output/n' num2str(N) 'TrueSimVesID' num2str(vesID) '_bgFlow' bgFlow '_speed' num2str(speed)];
else
fileName =  ['./output/n' num2str(N) 'TrueSimVesShape' initShape '_bgFlow' bgFlow '_speed' num2str(speed)];
end

% TIME STEPPING
for it = 2 : numel(timeTrue)
  disp('********************************************') 
  disp([num2str(it) 'th (/' num2str(numel(timeTrue)) ...
    ') time step, time: ' num2str(timeTrue(it))])
  
  tStart = tic;
  disp('---------------------------')
  if ~iSplit
    % SOLVE WITHOUT OPERATOR SPLITTING  
    disp('Solving without operator splitting...') 
    [XhistTrue(:,it),errALTrue(it)] = dnn.exactlySolve(XhistTrue(:,it-1),...
        vinf,dt,area0,len0,oc,op,kappa);
  else
    % SOLVE WITH OPERATOR SPLITTING  
    disp('Solving with operator splitting...')  
    % 1) TRANSLATION w/o NETWORK
    Xnew = dnn.translateVinf(vinf(XhistTrue(:,it-1)),dt,XhistTrue(:,it-1),op);   
    % 2) RELAXATION w/o NETWORK
    vesicle = capsules(Xnew,[],[],1,1,1);
    XhistTrue(:,it) = dnn.relaxExactSolve(vesicle,zeros(size(Xnew)),dt,Xnew,op);  
  end

  % AREA-LENGTH CORRECTION
  [Xnew,ifail] = oc.correctAreaAndLength2(XhistTrue(:,it),area0,len0);
  if ifail
    disp('Error in AL cannot be corrected!!!')
  end
  XhistTrue(:,it) = oc.alignCenterAngle(XhistTrue(:,it),Xnew);

  % check if shape is intersecting
  [xIntersect,~,~] = oc.selfintersect(XhistTrue(:,it));
  if ~isempty(xIntersect)
    disp('New vesicle shape is self-intersecting!!!')
  end
 
  if iequalDist
  % Equally distribute points in arc-length
  Xiter = XhistTrue(:,it);
  for iter = 1 : 5
    [Xiter,~,~] = oc.redistributeArcLength(Xiter);
  end
  % Fix misalignment in center and angle due to reparametrization
  XhistTrue(:,it) = oc.alignCenterAngle(XhistTrue(:,it),Xiter);
  elseif ireparam
  [Xiter,~] = oc.reparametrize(XhistTrue(:,it),-XhistTrue(:,it)-XhistTrue(:,it-1),6,20);
  XhistTrue(:,it) = oc.alignCenterAngle(XhistTrue(:,it),Xiter);
  end

  disp(['took ' num2str(toc(tStart)) ' seconds.'])
  disp('---------------------------')    
  
  % Compute error in area and length
  [~,area,len] = oc.geomProp(XhistTrue(:,it));
  errALTrue(it) = max(abs(area-area0)/area0,abs(len-len0)/len0);
  disp(['Error in area and length: ' num2str(errALTrue(it))])   
  disp('********************************************') 
  disp(' ')
  
  if rem(it,100) == 0
    save(fileName,'XhistTrue','dt','timeTrue','errALTrue','vesID')
    figure(1);clf;
    plot(XhistTrue(1:end/2,it), XhistTrue(end/2+1:end,it))
    axis equal
    pause(0.1)
  end
  
  if strcmp(bgFlow,'rotation')
    cx = mean(XhistTrue(1:end/2,it)); cy = mean(XhistTrue(end/2+1:end,it));
    if sqrt(cx^2+cy^2)<=3*vesRadius
      save(fileName,'XhistTrue','dt','timeTrue','errALTrue','vesID','it')
      disp('Can kill the simulation, vesicle approaches to origin')
    end
  end
end

% SAVE DATA 
save(fileName,'XhistTrue','dt','timeTrue','errALTrue','vesID')



