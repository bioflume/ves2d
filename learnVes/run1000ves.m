clear; clc;
addpath ../src/
addpath ../examples/
oc = curve;

chanWidth = 1.8;
speed = 1500;
timeHorizon = 0.5; 
dtGiven = 1e-4;

%-------------------------------------------------------------------------

% VESICLE:
% -------------------------------------------------------------------------

load FewFig8IC.mat
X = Xic;
[~,area0,len0] = oc.geomProp(X);
nv = numel(X(1,:));
N =  32; % num. points


disp('Single vesicle with background fluid, Exact Solve')

% PARAMETERS, TOOLS
%-------------------------------------------------------------------------
prams.repStrength = 0;
prams.bgFlow = 'parabolic'; % 'shear','tayGreen','relax','parabolic'
prams.speed = speed; % 500-3000 for shear, 70 for rotation, 100-400 for parabolic
prams.Th = timeHorizon;
prams.N = N; % num. points for true solve in DNN scheme
prams.Nfmm = N;
prams.nv = nv; %(24 for VF = 0.1, 47 for VF = 0.2) num. of vesicles
prams.fmm = false; % use FMM for ves2ves
prams.fmmDLP = false; % use FMM for ves2walls
prams.kappa = 1;
prams.dt = dtGiven; % time step size
prams.dtRelax = prams.dt;
prams.Nbd = 0;
prams.nvbd = 0;
prams.interpOrder = 1;
prams.chanWidth = chanWidth;
dnn = dnnToolsManyVesFree(X,prams);


% INITIALIZE MATRICES AND COUNTERS
% ------------------------------------------------------------------------
% folderName = '/work2/03353/gokberk/frontera/truePoisRuns/';
% fileName = ['./output/32modes_poisTrueRuns_dt' num2str(dt) '_speed' num2str(speed) '_width' num2str(chanWidth) '.bin'];
fileName = ['./output/test.bin'];

fid = fopen(fileName,'w');
output = [N;nv];
fwrite(fid,output,'double');

x = X(1:end/2,:); y = X(end/2+1:end,:);
output = [x(:); y(:)];
fwrite(fid,output,'double');
fclose(fid);

timeTrue = [0];
XhistTrue = X; sigStore = zeros(N,nv);

writeData(fileName,XhistTrue,sigStore,timeTrue(end),0,0);

% TIME STEPPING
it = 1;
cx = []; cy = [];
while timeTrue(end) < prams.Th
  disp('********************************************') 
  disp([num2str(it) 'th time step, time: ' num2str(timeTrue(it))])
  
  tStart = tic;
  disp('---------------------------')
  
  % SOLVE WITHOUT OPERATOR SPLITTING  
  disp('Solving without operator splitting...') 
  % SOLVE WITHOUT OPERATOR SPLITTING  

  [Xnew,sigStore] = dnn.DNNlikeExactSolve(XhistTrue,sigStore);

  
  
  % AREA-LENGTH CORRECTION
  [Xnew2,ifail] = oc.correctAreaAndLength2(Xnew,area0,len0);
  if ifail
    disp('Error in AL cannot be corrected!!!')
  else
    Xnew = oc.alignCenterAngle(Xnew, Xnew2);
  end
  
 
  
  XOrig = Xnew;
  for iter = 1 : 5
    Xnew = oc.redistributeArcLength(Xnew);
  end
  XhistTrue = oc.alignCenterAngle(XOrig,Xnew);

  it = it + 1;
  timeTrue(it) = timeTrue(it-1) + prams.dt;  


  figure(1);clf;
  plot(XhistTrue(1:end/2,:),XhistTrue(end/2+1:end,:),'r','linewidth',2)
  axis equal
  % ylim([-chanWidth, chanWidth])

  disp(['took ' num2str(toc(tStart)) ' seconds.'])
  disp('---------------------------')    
  
  % Compute error in area and length
  [~,area,len] = oc.geomProp(XhistTrue);
  errALTrue = max(abs(area-area0)/area0,abs(len-len0)/len0);
  disp(['Error in area and length: ' num2str(errALTrue)])   
  disp('********************************************') 
  disp(' ')
  
  

  if rem(it,10) == 0
    writeData(fileName,XhistTrue,sigStore,timeTrue(end),0,0);  
 
  end
  
end % while


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeData(filename,X,sigma,time,ncountNN,ncountExact)
x = X(1:end/2,:);
y = X(end/2+1:end,:);
output = [time;ncountNN;ncountExact;x(:);y(:);sigma(:)];

fid = fopen(filename,'a');
fwrite(fid,output,'double');
fclose(fid);


end




